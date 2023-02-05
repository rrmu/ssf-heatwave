# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 23:24:37 2022

@author: Michael
"""
import netCDF4 as nc
import numpy as np

def toNC(file,exportname,unit):
    ds = nc.Dataset(exportname, 'w', format='NETCDF4')
    
    ds.createDimension('time', file[:,1,1].size)
    ds.createDimension('lat', 360)
    ds.createDimension('lon', 720)
    
    times = ds.createVariable('time', 'i', ('time',))
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))
    value = ds.createVariable('value', 'f4', ('time', 'lat', 'lon',))
    value.units = unit
    
    times[:] = np.arange(1,file[:,1,1].size)
    lats[:] = np.arange(89.75, -89.75, 0.5)
    lons[:] = np.arange(0.25, 359.75, 0.5)
    value[:,:,:]=file[:,:,:]
    
    ds.close()