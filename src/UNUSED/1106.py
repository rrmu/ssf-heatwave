# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 19:16:45 2022

@author: Michael
"""

import netCDF4 as nc
import numpy as np
from progress.bar import Bar

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
        # file是个字典，其中"tmax"键值对是一个三维数组，形状至少为（file.variables['time'].size，360，720）
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