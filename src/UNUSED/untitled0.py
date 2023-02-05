# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:59:58 2022

@author: Michael
"""
class Dataset:
    def getheatcase(path, percentileval):
        import numpy as np
        import netCDF4 as nc
        
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