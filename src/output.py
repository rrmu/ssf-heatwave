# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 16:42:59 2022

Graph data

@author: Michael
"""
import os
import rioxarray  # NOTE: this package is called with alias xarray.Dataset.rio

import numpy as np
import xarray as xr
import geopandas as gpd
import netCDF4 as nc
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.ticker import MultipleLocator
from scipy.io import netcdf_file


class cdf:
    def percentile2cdf(origfile, p, name):
        """
        Used to export 2D numpy array over a single time frame.

        Parameters
        ----------
        origfile : 3D/2D NetCDF4 File
            Where your percentile/eventdata is from.
        p : 2D NumPy Array
            Percentile Array.
        name : New Filename

        Exports
        -------
        File.

        """

        file = nc.Dataset(origfile)
        time = np.array([1000])
        lat_arr = file.variables['lat'][:].data
        lon_arr = file.variables['lon'][:].data
        filename = netcdf_file(name+'.nc', 'w')
        ntime, nlat, nlon = len(time), len(lat_arr), len(lon_arr)

        # Dimensions
        filename.createDimension('time', ntime)
        filename.createDimension('lat', nlat)
        filename.createDimension('lon', nlon)

        # Variables
        time = filename.createVariable('time', 'i', ('time',))
        lat = filename.createVariable('lat', 'f4', ('lat',))
        lon = filename.createVariable('lon', 'f4', ('lon',))
        tmax = filename.createVariable('TmaxThreshold',
                                       'f4',
                                       ('time', 'lat', 'lon',)
                                       )

        # Attributes
        time.units = 'across time'
        lat.units = 'degrees north'
        lon.units = 'degrees east'
        tmax.units = 'degC'

        # Populate the variables with data
        time[:] = time.data
        lat[:] = lat_arr
        lon[:] = lon_arr
        tmax[:] = p[:]

        filename.close()

    def event2cdf(reffile, years_tuple, mod='exp', vname='var'):
        """
        Parameters
        ----------
        reffile : Path to 3D NetCDF4 File
            Same Lat Lon accuracy as where your data is from.
            This file is used only for lat lon reference.
         years_tuple : Tuple(numyear,3D NumPy Array)
             Array Tuple.
         mod : String
             Export filename modifier
        Exports
        -------
        Netcdf File.

        """

        file = nc.Dataset(reffile)

        def __stack(tupset):
            """    
            Tuple List MUST NOT HAVE 'Time' DIMENSION
            """
            arraylist = []
            for time in tupset:
                (year, array) = time
                arraylist.append(np.ma.masked_array(array,
                                                    mask=file.variables['lat'][1, :].mask)
                                 )

            stack = np.ma.stack(arraylist, axis=0)
            return stack

        t = []
        for year in years_tuple:
            (numyear, array) = year
            t.append(numyear)

        time = np.array(t)
        lat_arr = file.variables['lat'][:].data
        lon_arr = file.variables['lon'][:].data
        filename = netcdf_file(mod
                               + str(years_tuple[0][0])
                               + str(years_tuple[-1][0])
                               + '.nc',
                               'w')
        ntime, nlat, nlon = len(time), len(lat_arr), len(lon_arr)

        # Dimensions
        filename.createDimension('time', ntime)
        filename.createDimension('lat', nlat)
        filename.createDimension('lon', nlon)

        # Variables
        time = filename.createVariable('time', 'i', ('time',))
        lat = filename.createVariable('lat', 'f4', ('lat',))
        lon = filename.createVariable('lon', 'f4', ('lon',))
        tmax = filename.createVariable(vname, 'f4', ('time', 'lat', 'lon',))

        # Attributes
        time.units = 'year'
        lat.units = 'degrees north'
        lon.units = 'degrees east'
        tmax.units = 'degC'

        # Populate the variables with data
        time[:] = time
        lat[:] = lat_arr
        lon[:] = lon_arr
        tmax[:, :, :] = __stack(years_tuple)

        filename.close()


class Continental:
    def __init__(self,
                 cdf,
                 shp,
                 name,
                 shp_is_file=False,
                 lon_shift=False,
                 maskref='/Users/maeko/Documents/ISEF2023_EAEV_011/Temperature/1/tmax.1980.nc',
                 dims=['time', 'lat', 'lon']
                 ):
        """


        Parameters
        ----------
        ds : NetCDF4 file PATH
            Originial dataset.
        shapefile : shapefile PATH
            Path of desired shapefile.
        lon_shift : BOOL, optional
            Shift longitude to (-180,180). The default is False.
        dimslist : List of str, optional
            List the name of dimensions. The default is ['time', 'lat', 'lon'].

        Raises
        ------
        ValueError
            An exception that occurs when function receives an argument of the correct data type but an inappropriate value.
        TypeError
            An exception that occurs when the data type of an object in an operation is inappropriate..

        Returns
        -------
        Xarray.Dataset
            Geographically clipped xarray dataset of the original data.

        """

        def clip_geo(ds, shape, dimslist,conv_back=False,dropnan=True):
            if shp_is_file:
                sf = gpd.read_file(shape)
            else:
                sf = gpd.GeoDataFrame(
                    {'col1': ['name1'], 'geometry': [shape]}, crs="EPSG:3857")
            """
            def __findcoords(array):
                __coorddict = {}
                for dim in dimslist:
                    __coorddict[dim] = range(__da[dim][0], __da[dim][-1])
            """
            # Read data from ds as xarray

            print('Getting data...')
            try:
                # if type(ds) is str:
                reffile = nc.Dataset(maskref)

                msk = reffile.variables['tmax'][0:43, :].mask

                origds = xr.open_dataset(ds)
                a = origds.HeatWaveOccurence.to_masked_array()

                marray = np.ma.masked_where(msk, a)

                dimstup = tuple(dimslist)
                __da = xr.Dataset(data_vars={'HeatWaveOccurence': (dimstup,
                                                                   marray.filled(
                                                                       fill_value=np.nan)
                                                                   )
                                             },
                                  coords={'time': origds['time'],
                                          'lat': origds['lat'],
                                          'lon': origds['lon'],
                                          }
                                  )

            except:
                raise TypeError('Data type '
                                + "'"
                                + str(type(ds).__name__)
                                + "'"
                                + ' is not endorsed.')

            # Clip dataset
            else:

                print('Clipping data...')
                __da.rio.write_crs("EPSG:3857", inplace=True)

                # longitude shift

                if lon_shift:
                    __da.coords['lon'] = (__da.coords['lon'] + 180) % 360 - 180
                    __da = __da.sortby(__da.lon)
                
                
                # clip geometry
                self.original = __da
                __da.rio.set_spatial_dims(x_dim='lon', y_dim='lat')
                clipped=__da.rio.clip(sf.geometry, __da.rio.crs, drop=dropnan)
                if conv_back:
                    clipped.coords['lon'] = (clipped.coords['lon'] +180)
                    clipped = clipped.sortby(clipped.lon)
                    
                    return clipped
                else: return clipped

        self.data = clip_geo(cdf, shp, dims)  # xarray dataset
        self.raw = clip_geo(cdf,shp,dims,conv_back=True,dropnan=False)
        self.name = name
        self.ma = np.ma.masked_invalid(self.raw.HeatWaveOccurence.to_numpy()).mask
        self.avr = self.data.HeatWaveOccurence.mean(
            dim={'lat', 'lon'}, skipna=True)
        self.med = self.data.HeatWaveOccurence.median(
            dim={'lat', 'lon'}, skipna=True)
        self.sd = self.data.HeatWaveOccurence.std(
            dim={'lat', 'lon'}, skipna=True)

    def save(self):
        parent = Path(os.getcwd()).parent
        __filename = os.path.join(parent, 'output/regional/')
        np.savetxt(__filename+self.name+"stats.csv", [self.avr,
                                                      self.med,
                                                      self.sd],
                   delimiter=",")

    def plot(self, t):
        f, ax = plt.subplots(1, 1)
        ax.pcolormesh(self.data.HeatWaveOccurence[t, :])
        ax.invert_yaxis()
        ax.axis('equal')
        ax.axis('off')


class YearStats:
    def __init__(self, UserTupList):
        self.pure = UserTupList
        self.avr = []
        self.med = []
        self.label = []
        self.sd = []

        for (year, marray) in UserTupList:
            self.avr.append(np.average(marray.compressed()))
            self.med.append(np.median(marray.compressed()))
            self.sd.append(np.nanstd(marray.compressed()))
            self.label.append(year)

    def box(self, title=None, *args):
        arr1d_li = []
        for (year, marray) in self.pure:
            arr1d_li.append(marray.compressed())

        plt.figure(figsize=(30, 9))

        plt.boxplot(arr1d_li,
                    labels=self.label,
                    showfliers=False,
                    showmeans=True,
                    meanline=True)
        if title == None:
            plt.title('Plot')
        else:
            plt.title(str(title))

        plt.show()

    def save(self):
        parent = Path(os.getcwd()).parent
        __filename = os.path.join(parent, 'output/global/')
        np.savetxt(__filename+"global_stats.csv", [self.label,
                                                   self.avr,
                                                   self.med,
                                                   self.sd],
                   delimiter=",")

    def line(self, title=None, export_data=False, *args):

        fig, ax = plt.subplots()
        ax.plot(self.label, self.avr)
        ax.plot(self.label, self.med)

        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter('{x:.0f}')

        # For the minor ticks, use no labels; default NullFormatter.
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))

        fig.set_dpi(200)
        plt.xlabel('Year')
        plt.ylabel('Heat Wave Occurence')

        if title == None:
            plt.title('Plot')
        else:
            plt.title(str(title))

        plt.plot()
        plt.show()

    def lat_avr(self, startcoord=-89.75, endcoord=89.75):
        arr1d_li = []
        fixed = []

        __label = []
        i = 0

        step = (endcoord-startcoord)/(len(self.pure[1][1][:, 1])-1)
        while(i < len(self.pure[1][1][:, 1])):
            __label.append(startcoord+step*i)
            i += 1

        for tup in self.pure:
            fixed.append(np.ma.fix_invalid(tup[1]))

        for l in fixed:
            arr1d_li.append(np.average(l.data, axis=0))

        plt.figure(figsize=(16, 9))

        plt.plot(arr1d_li)
        plt.xlabel('Latitude')
        plt.ylabel('Mean Heat Wave Occurence')
        plt.yticks(range(1, 10))

        plt.title('Plot')

        plt.plot()

        plt.show()
