# -*- coding: utf-8 -*-
"""
Calculate multidimentional time series

@author: Michael
"""

import netCDF4 as nc
import numpy as np
import xarray as xr
import scipy.stats as sps
from progress.bar import Bar

def percentile(Filename, percentile):
        """
        这个是计算Filename.tmax[lat, lon, time]每一个[lat, lon]上的percentile值
        返回一个[lat, lon]大小的矩阵
        """
        ds = xr.open_dataset(Filename)
        tmax=ds['tmax']
    
        mu = tmax.mean(dim=['time'],skipna=True)
        sigma = tmax.std(dim=['time'],skipna=True)
    
        dist = sps.norm(loc=mu, scale=sigma)
        p = dist.ppf(percentile)
        
        return p

def getheatcase(path, percentileval):
        """
    

    Parameters
    ----------
    path : String
        Must point to a 3D netCDF file following var[time,lat,lon].
    percentileval : NumPy Array
        Must be a 2D array.

    Returns
    -------
    [time,lat,lon]
    NumPy Array with no Attributes
        Consisting of 3 dimensions and bool values.

    """
        # file是个字典，其中"tmax"键值对是一个三维数组
        file=nc.Dataset(path)
        isheatday=np.zeros((file.variables['time'].size,360, 720))
    
        # i至少为360，j至少为720
        i=file.variables['lat'].size
        v=file.variables['lon'].size
        t = file.variables['time'].size
    
        data = file.variables['tmax'][:t,:i,:v]
    
        # index为提取所有not nan值的索引
        index = ~np.isnan(data[0])
    
        # 得到在percentileval的需要需要比较的对应值
        cmp_values_1d = percentileval[:i,:v][index]
    
        # 得到所有需要与percentileval比较的对应的数据
        cmp_values_2d = data[:,index]
    
        # 将2组数据进行比较，得到比较结果
        cmp_ret_2d = (cmp_values_2d > cmp_values_1d)
    
        # 进行赋值并返回
        isheatday[:,:i,:v][:,index] = cmp_ret_2d
        
        return np.array(isheatday,dtype=np.bool)
    
def lessdayremover(origset,minday):
    """
    

Parameters
----------
origset : NumPy Array
    Must point to a 3D Numpy Array.
minday : int
    Minimum accepted consecutive days.

Returns
-------
[lat,lon,time]
NumPy Array with no Attributes
    Consisting of 3 dimensions and bool values.

"""
    dataset=origset
    time=0
    thrownvalue=0
    totaltrue=0
    with Bar('Processing...', max=360) as bar:
        for lat in range(dataset[1,:,1].size):
                for lon in range(dataset[1,1,:].size):

                        counter=0
                        for time in range(dataset[:,1,1].size):
                            if dataset[time,lat,lon]:
                                totaltrue+=1
                                counter+=1
                            #If False
                            elif counter<minday:
                                for i in range(1,counter+1):
                                    dataset[time-i,lat,lon]=False
                                    totaltrue-=1
                                    thrownvalue+=1
                                counter=0
                            else: 
                                counter=0
                bar.next()                        
        print("")
        print("Values Thrown: ", thrownvalue, " out of ", dataset.size)
        print("Current True Values: ", totaltrue, " out of ", dataset.size)
        
        return dataset

def daycounter(dataset):
    """
    

    Parameters
    ----------
    file : 3D BOOL NumPy Array

    Returns
    -------
    [lat,lon]
    Integer NumPy Array with no Attributes

    """
    numday=np.zeros(360, 720)
    with Bar('Processing...', max=360) as bar:
        for lat in range(dataset[1,:,1].size):
                for lon in range(dataset[1,1,:].size):
                    for time in range(dataset[:,1,1].size):
                        if dataset[time,lat,lon]:
                            numday[lat,lon]+=1
                bar.next()
                    
                    

def convdata(origfile,p,name):
        from scipy.io.netcdf import netcdf_file
        import netCDF4 as nc
        
        file=nc.Dataset(origfile)
        time=np.array([1000])
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
        tmax = filename.createVariable('TmaxThreshold', 'f4', ('time','lat', 'lon',))
        
        # Attributes
        time.units = 'across time'
        lat.units = 'degrees north'
        lon.units = 'degrees east'
        tmax.units = 'degC'
        
        # Populate the variables with data
        time[:]=time.data
        lat[:] = lat_arr
        lon[:] = lon_arr
        tmax[:,:,:] = p[:,:].data
        
        filename.close()