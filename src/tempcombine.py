# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 12:52:12 2022

@author: Michael
"""

import xarray

ds1 = xarray.open_mfdataset('D:/Temperature/refbank/tmax.*.nc',combine = 'nested', concat_dim="time")
ds1.to_netcdf('D:/Temperature/tmax_climatology_8110.nc')