#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 21:04:15 2022

@author: maeko
"""
import numpy as np
import netCDF4 as nc
import xarray as xr

reffile=nc.Dataset('/Users/maeko/Documents/ISEF2023_EAEV_011/Temperature/4/tmax.2021.nc')

msk=reffile.variables['tmax'][0:43,:].mask

origds=xr.open_dataset('hwcase19802022.nc')
a=origds.HeatWaveOccurence.to_masked_array()

marray=np.ma.masked_where(msk,a)

dimstup=('time','lat','lon')
newds = xr.Dataset(data_vars={'HeatWaveOccurence': (dimstup,
                                          marray.filled(fill_value=np.nan)
                                          )
                    },
                   coords={'time': origds['time'],
                           'lat': origds['lat'],
                           'lon': origds['lon'], 
                           }
                   )
