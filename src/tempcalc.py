# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 21:05:52 2022

@author: Michael

EXECUTION PROGRAM FOR HEATWAVE ANALYSIS
"""
import os
import concurrent.futures
import tqdm
import output

import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as sps
import netCDF4 as nc

from pathlib import Path
from shapely.geometry import Polygon

#------Heatwave Algorithm------#


def get_percentile(Filename, percentile):
    """


    Parameters
    ----------
    Filename : 3D NetCDF Dataset
        Data source containing a ['time'] dimension.
    percentile : Float, 0<value<1
        Desired percentile (left to right).

    Returns
    -------
    p : 2D NumPy array
        90th percentile of Filename values across time.

    """

    ds = xr.open_dataset(Filename)
    tmax = ds['tmax']

    # Getting necessary descriptives
    mu = tmax.mean(dim=['time'], skipna=True)
    sigma = tmax.std(dim=['time'], skipna=True)

    # Generate distribution
    dist = sps.norm(loc=mu, scale=sigma)
    p = dist.ppf(percentile)

    return p


def caseof(path, percentileval, minday, casecount=True, display=False):
    """


    Parameters
    ----------
    path : String
        Must point to a 3D netCDF file following var[time,lat,lon].
    percentileval : NumPy Array
        Must be a 2D array.

    Returns
    -------
    casecount=False
    (int year,arraydims[time,lat,lon])
    NumPy Array with no Attributes
        Consisting of 3 dimensions and bool values.
    casecount=True
    (int year,arraydims[lat,lon])

    """
    def __lessdayremover(origset, minday):

        time = 0
        thrownvalue = 0
        totaltrue = 0

        dataset = origset
        eventcount = np.zeros((360, 720))

        for lat in range(dataset[1, :, 1].size):
            for lon in range(dataset[1, 1, :].size):

                counter = 0
                for time in range(dataset[:, 1, 1].size):
                    """
                        Time pointer on each lat lon.
                        If is True, add 1 to counter;
                        If is False, check if counter<minday,
                        If so, delete all TRUEs counted.
                        Each heat case end with counter>3, then counter=0
                        """
                    if dataset[time, lat, lon]:
                        totaltrue += 1
                        counter += 1
                    # False counter<minday
                    elif dataset[time, lat, lon] == False:
                        if counter < minday:
                            for i in range(1, counter + 1):
                                dataset[time - i, lat, lon] = False
                                totaltrue -= 1
                                thrownvalue += 1
                            counter = 0
                        # False counter>=minday, eventcount[lat lon]++
                        else:
                            eventcount[lat, lon] += 1
                            counter = 0
        if display:
            print("")
            print(path, ":")
            print("Values Thrown: ", thrownvalue, " out of ", dataset.size)
            print("Current True Values: ", totaltrue, " out of ", dataset.size)

        print()
        print("Array ready")

        res_mask = np.ma.getmask(data)

        if casecount:
            return np.ma.masked_array(eventcount, mask=res_mask[0])
        else:
            return np.ma.masked_array(dataset, mask=res_mask)

    print("Acquiring Heat Day")
    
    file = nc.Dataset(path)
    isheatday = np.empty((file.variables['time'].size, 360, 720))
    isheatday[:] = np.nan

    # min i=360，min j=720
    i = file.variables['lat'].size
    v = file.variables['lon'].size
    t = file.variables['time'].size

    data = file.variables['tmax'][:t, :i, :v]

    # indexing all non NaN values
    index = ~np.isnan(data[0])

    # get corresponding values on percentileval
    cmp_values_1d = percentileval[:i, :v][index]

    # 得到所有需要与percentileval比较的对应的数据
    cmp_values_2d = data[:, index]

    # 将2组数据进行比较，得到比较结果
    cmp_ret_2d = (cmp_values_2d > cmp_values_1d)

    # 进行赋值并返回
    isheatday[:, :i, :v][:, index] = cmp_ret_2d

    print("Analyzing Data")

    return (int(os.path.basename(path).split(".")[1]),
            __lessdayremover((np.array(isheatday, dtype=np.bool)), minday))


#------Multithreading------#


def get_case(directory, p, mindays, isevent=False, display=False):
    """


    Parameters
    ----------
    directory : LIST
        List of COMPLETE directories.
    p : NUMPY ARRAY
        Percentile.
    mindays : INT
        Minimum day to be considered heat wave.

    Returns
    -------
    yearboollist : LIST of TUPLES
        List of tuples (Year, NumPy bool arrays).

    """
    yearboollist = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=(16)) as e:
        futures = [e.submit(caseof, dataset, p, mindays, isevent, display)
                   for dataset in directory]
        for group in concurrent.futures.as_completed(futures):
            yearboollist.append(group.result())

    if isevent:
        print('Event count finished for list of', len(yearboollist), 'arrays.')
    else:
        print('Computation process finished for list of',
              len(yearboollist), 'arrays.')
    return sorted(yearboollist)


#------Supplementary Methods------#


def lats(li, axis=0):
    #res_mask = np.ma.getmask(ds[0])
    tup = []
    for tuplist in li:
        for (year, array) in tuplist:
            tup.append(np.average(array, axis=1))
    __lats = np.ma.stack(tup, axis=0)
    return __lats.filled(np.nan)


def to_tupset(netcdf_path,
              start_year,
              mask_reference=None,
              var_name='HeatWaveOccurence'):

    __ds = nc.Dataset(netcdf_path)

    tuplist = []
    if mask_reference != None:
        reffile = mask_reference
    else:
        reffile = netcdf_path

    for year in range(len(__ds.variables['time'][:].data)):
        tuplist.append((start_year+year,
                        __ds.variables[var_name][year, :]
                        )
                       )
        print(tuplist[year][0])

    return masking(tuplist, reffile, var=var_name)


def masking(origlist, geo_nan_file_path, var='tmax'):

    origds = nc.Dataset(geo_nan_file_path)
    newds = origds.variables[var][:]

    masked = []
    res_mask = np.ma.getmask(newds[0])
    for i in range(len(origlist)):

        if type(origlist[i]) is tuple:
            masked.append((origlist[i][0],
                           np.ma.masked_array(
                               origlist[i][1].data, mask=res_mask)
                           ))

        elif type(origlist[i]) is np.ndarray:
            masked.append(np.ma.masked_array(
                origlist[i].data, mask=res_mask)
            )

        else:
            raise TypeError

    return masked

#------Execution------#


if __name__ == '__main__':
    parent = Path(os.getcwd()).parent
    grandparent = Path(parent).parent
# OBSERVED WORLD DATA
    reffile = os.path.join(
        grandparent, 'Temperature/tmax_climatology_8110.nc').replace("\\", "/")

    print("===PREPARING DATA===")

    # GET DIRECTORY LIST (OBSERVED DATA)
    def completedir(dirname):
        newdir = []
        for file in os.listdir(dirname):
            newdir.append(os.path.join(dirname, file).replace("\\", "/"))
        return newdir

    dir1 = completedir(os.path.join(grandparent, 'Temperature/1'))
    dir2 = completedir(os.path.join(grandparent, 'Temperature/2'))
    dir3 = completedir(os.path.join(grandparent, 'Temperature/3'))
    dir4 = completedir(os.path.join(grandparent, 'Temperature/4'))

    # GET PERCENTILE
    percentile = get_percentile(reffile, .9)

    print("===ANALYZING HEATWAVE CASE===")

    # GET HEATWAVE EVENT
    yearevent8089 = get_case(dir1, percentile, 3, isevent=True)
    yearevent9099 = get_case(dir2, percentile, 3, isevent=True)
    yearevent0009 = get_case(dir3, percentile, 3, isevent=True)
    yearevent1022 = get_case(dir4, percentile, 3, isevent=True)

    # ORGANIZE RESULTS
    yearevent = []
    for tup in yearevent8089:
        yearevent.append(tup)

    for tup in yearevent9099:
        yearevent.append(tup)

    for tup in yearevent0009:
        yearevent.append(tup)

    for tup in yearevent1022:
        yearevent.append(tup)

    # EXPORT WORLD DATA
    output.percentile2cdf(reffile, percentile, 'percentile90_climatology')
    output.event2cdf(reffile, yearevent, 'hwcase', 'HeatWaveOccurence')
