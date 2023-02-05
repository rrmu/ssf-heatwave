#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 16:36:26 2022

@author: maeko
"""
import os
import output
import pandas as pd
import numpy as np
import xarray as xr
import netCDF4 as nc
import geopandas as gpd
import matplotlib.pyplot as plt
import tqdm
from pathlib import Path
from shapely.geometry import Polygon

parent = Path(os.getcwd()).parent
grandparent = Path(parent).parent

# REGIONAL HEATWAVE MEASURES
# CONTINENTAL SHAPEFILE PATHS
eu_ru = Polygon(((-10.5, 36),
                (180, 36),
                 (180, 87.75),
                 (-10.5, 87.75)
                 ))

n_a = Polygon(((-170, 32.5),
               (-11.5, 32.5),
               (-11.5, 87.75),
               (-170, 87.75)
               ))

la_ca = Polygon(((-117.25, -56),
                (-34.75, -56),
                 (-34.75, 32.5),
                 (-117.25, 32.5)
                 ))

a_p = Polygon(((35, -50),
               (180, -50),
               (180, 53.5),
               (35, 53.5),
               ))

afr = os.path.join(
    parent,
    'shapefile/africa/icpac-afr_g/').replace("\\", "/")

# GET CONTINENTAL HEATWAVE EVENT COUNTS
d_africa = output.Continental(
    'hwcase19802022.nc', afr, 'afr', True, True)
d_europe = output.Continental(
    'hwcase19802022.nc', eu_ru, 'eu_ru', False, True)
d_n_america = output.Continental(
    'hwcase19802022.nc', n_a, 'n_a', False, True)
d_s_america = output.Continental(
    'hwcase19802022.nc', la_ca, 'la_ca', False, True)
d_ap = output.Continental(
    'hwcase19802022.nc', a_p, 'a_p', False, True)

year=41 # Taken Year of 2021
"""
d_africa.plot(year)
d_europe.plot(year)
d_n_america.plot(year)
d_s_america.plot(year)
d_ap.plot(year)

plt.show()

d_africa.save()
d_europe.save()
d_n_america.save()
d_s_america.save()
d_ap.save()
"""
# Construct Correlation Map Basemap

# Create maskedarray
lpi_size=(49,360,720)

# Global
file_glob = pd.read_csv(
    os.path.join(grandparent,'LPI/Global.csv')
    )

glob=[]
for thing in file_glob['LPI_final']:
    glob.append(thing)

tmax_mask_path=os.path.join(grandparent, 'Temperature/4/tmax.2021.nc')
with nc.Dataset(tmax_mask_path) as reffile:
    msk=reffile.variables['tmax'][0,:].mask
    

# Read CSVs
# Africa
file_afr=pd.read_csv(
    os.path.join(grandparent,'LPI/Africa.csv')
    )

africa=[]
for thing in file_afr['LPI_final']:
    africa.append(thing)

# North America
file_n_a=pd.read_csv(
    os.path.join(grandparent,'LPI/North America.csv')
    )

n_america=[]
for thing in file_n_a['LPI_final']:
    n_america.append(thing)

# Latin America and Caribbeans
file_laca=pd.read_csv(
    os.path.join(grandparent,'LPI/Latin America and Caribbean.csv')
    )

s_america=[]
for thing in file_laca['LPI_final']:
    s_america.append(thing)

# Europe and Central Asia
file_euru=pd.read_csv(
    os.path.join(grandparent,'LPI/Europe and Central Asia.csv')
    )

eu_ru=[]
for thing in file_euru['LPI_final']:
    eu_ru.append(thing)
    
# Asia Pacific
file_ap=pd.read_csv(
    os.path.join(grandparent,'LPI/Asia and the Pacific.csv')
    )

a_p=[]
for thing in file_ap['LPI_final']:
    a_p.append(thing)
    


order = [
         (d_ap.ma[1],a_p), 
         (d_n_america.ma[1],n_america), 
         (d_s_america.ma[1],s_america), 
         (d_europe.ma[1],eu_ru), 
         (d_africa.ma[1],africa)]

# order[n][0] -> bool numpy array of masks
# order[n][1] -> pandas 1d dataframe (time series)
# Each mask has same (time, lat, lon) as in "order"


#correlate over 'order'
correl=np.zeros((360,720))
ds=nc.Dataset('hwcase19802022.nc')
heatwave=ds.variables['HeatWaveOccurence'][:].data

for area in order:
    (mask,t_list)=area
    t_series=pd.Series(t_list)
    for i in tqdm.tqdm(range(100), desc="Calculating..."):
        for lat in range(len(heatwave[0,:,0])):
            for lon in range(len(heatwave[0,0,:])):
                if mask[lat,lon]:
                        __dfgeo=pd.Series(heatwave[:38,lat,lon])
                        correl[lat,lon]=__dfgeo.corr(t_series)
    
plt.pcolormesh(correl)
plt.colorbar()
plt.gca().invert_yaxis()

lon_p=np.arange(0.25,180.25,0.5)
lon_n=np.arange(-179.75,0.25,0.5)
lon_coef=np.append(lon_p,lon_n)
lat_c=reffile.variables['lat'][:].data-90
lat_c=np.flip(lat_c)

# Export correlation map
ds=xr.Dataset({'PearsonCoef':(["lat","lon"],correl)},
              coords={
                      "lon": (["lon"], lon_coef),
                      "lat": (["lat"], lat_c),
                      "time":(['time'], [1]),
                      })

ds = ds.sortby(ds.lon)