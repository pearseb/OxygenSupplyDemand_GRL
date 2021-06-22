# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 11:44:53 2017

@author: pearseb
"""

#%% imports

from __future__ import unicode_literals

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import GridSpec
import netCDF4 as nc
import cmocean.cm as cmo
from scipy.optimize import curve_fit
from matplotlib.animation import ArtistAnimation
import seaborn as sb
sb.set(style='ticks')
from tqdm import tqdm
import pandas as pd


#%% get GLODAPv2 (2019) d13c data of DIC

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data')
glodap = pd.read_csv('GLODAPv2.2020_Merged_Master_File.csv', usecols=(3,4,8,9,14,15,17,26))

print(glodap.columns)
print(np.min(glodap['year']), np.max(glodap['year']))
print(np.min(glodap['month']), np.max(glodap['month']))
print(np.min(glodap['longitude']), np.max(glodap['longitude']))
print(np.min(glodap['oxygen']), np.max(glodap['oxygen']))

### remove rows with nan values
glodap = glodap.replace(-9999.0, np.nan)
glodap = glodap.dropna()


### make all longitudes positive
glodap['longitude'][glodap['longitude']<0] += 360.0

print("Number of GLODAPv2 (2019) oxygen measurements =", len(glodap['oxygen']))



#%% load model grid

os.chdir('C:\\Users\\pearseb\\Dropbox\\ARISE\\Modelling\\assessment')
data = nc.Dataset('ETOPO_ORCA2_OFF_PISCESiso_spinup_1y_ptrc_Y4200.nc', 'r')
no3 = data.variables['NO3'][0,...]
lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]
dep = data.variables['deptht'][...]

# fix longitudes as all positive
lon[lon>360.0] -= 360.0


### get bounds of coordinate grid
dep_bnds = data.variables['deptht_bnds'][...]
lon_bnds = np.zeros((len(lon),2))
lat_bnds = np.zeros((len(lat),2))
lon_bnds[:,0] = lon[:]-0.5; lon_bnds[:,1] = lon[:]+0.5
lat_bnds[:,0] = lat[:]-0.5; lat_bnds[:,1] = lat[:]+0.5

data.close()


#%% average data at each grid point of the model and do so for every year

years = np.arange(np.min(glodap['year']), np.max(glodap['year'])+1, 1)

oxy_grid = np.zeros((len(years),12,31,180,360))
count = np.zeros(np.shape(oxy_grid))

# for every year
for yr,year in enumerate(years):
    # for each month 
    for mn,mnth in enumerate(np.arange(1,13,1)):
        # select only data in month
        tmp = glodap.loc[glodap['year'] == year]
        tmp = tmp.loc[tmp['month'] == mnth]
        if len(tmp['oxygen']) == 0:
            oxy_grid[yr,mn,:,:,:] = np.nan
            continue
        
        print("YEAR", year, "MONTH", mnth, "DATA", len(tmp['oxygen']))
        tmp1 = tmp
        
        # for each grid cell
        for i in tqdm(np.arange(len(lon)), desc='Longitudes'):
            #print(i,lon_bnds[i,0], lon_bnds[i,1])
            tmp2 = tmp1.loc[ (tmp1['longitude'] >= lon_bnds[i,0]) & (tmp1['longitude'] < lon_bnds[i,1]) ]
            if len(tmp2['oxygen']) > 0:
                for j in np.arange(len(lat)):
                    #print(j,lat_bnds[j,0], lat_bnds[j,1])
                    tmp3 = tmp2.loc[ (tmp2['latitude'] >= lat_bnds[j,0]) & (tmp2['latitude'] < lat_bnds[j,1]) ]
                    if len(tmp3['oxygen']) > 0:
                        for k in np.arange(len(dep)):
                            #print(k,dep_bnds[k,0], dep_bnds[k,1])
                            tmp4 = tmp3.loc[ (tmp3['depth'] >= dep_bnds[k,0]) & (tmp3['depth'] < dep_bnds[k,1]) ]
                            if len(tmp4['oxygen']) > 0:
                                print("Number of observations at year", year, "and month", mnth, "after location selection =", len(tmp4['oxygen']))
                                oxy_grid[yr,mn,k,j,i] = np.mean(tmp4['oxygen'])
                                count[yr,mn,k,j,i] = len(tmp4['oxygen'])
                            else:
                                oxy_grid[yr,mn,k,j,i] = np.nan
                    else:    
                        oxy_grid[yr,mn,:,j,i] = np.nan
            else:
                oxy_grid[yr,mn,:,:,i] = np.nan
                
print(np.nanmin(oxy_grid), np.nanmax(oxy_grid))


for yr,year in enumerate(years):
    for mn,mnth in enumerate(np.arange(1,13,1)):
        print("Year =", year, "Month =", mnth, np.nanmin(oxy_grid[yr,mn,:,:,:]), np.nanmax(oxy_grid[yr,mn,:,:,:]))


#%% save as numpy array

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\oxygen')
np.savez('GLODAPv2.2020_oxygen_(yr,mn,dep,lat,lon)_1x1degrees.npz', oxy_grid)


#%% save as netcdf

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\oxygen')

os.remove('GLODAPv2.2020_oxygen_(yr,mn,dep,lat,lon)_1x1degrees.nc')
data = nc.Dataset('GLODAPv2.2020_oxygen_(yr,mn,dep,lat,lon)_1x1degrees.nc', 'w', format='NETCDF4_CLASSIC')

yrd = data.createDimension('t', len(years))
mnd = data.createDimension('m', 12)
xd = data.createDimension('x', 360)
yd = data.createDimension('y', 180)
zd = data.createDimension('z', 31)

yearv = data.createVariable('year', np.float64, ('t',))
mnthv = data.createVariable('month', np.float64, ('m',))
lonv = data.createVariable('lon', np.float64, ('x',))
latv = data.createVariable('lat', np.float64, ('y',))
depv = data.createVariable('dep', np.float64, ('z',))
oxyv = data.createVariable('oxygen', np.float64, ('t', 'm', 'z', 'y', 'x'))
countv = data.createVariable('count', np.float64, ('t', 'm', 'z', 'y', 'x'))

data.description = 'regridded oxygen values from GLODAPv2 2020'
data.history = "Created by Pearse J. Buchanan on 11th Dec 2020"
data.source = "GLODAPv2 (2020) downloaded from https://www.nodc.noaa.gov/ocads/oceans/GLODAPv2_2020/"

yearv.units = "years since 0001-01-01 00:00:00"
mnthv.units = "month"
lonv.units = "degrees_east"
latv.units = "degrees_north"
depv.units = "metres"
oxyv.units = "umol/kg"
countv.units = "number of measurements"

yearv.standard_name = "year"
mnthv.standard_name = "month"
lonv.standard_name = "longitude"
latv.standard_name = "latitude"
depv.standard_name = "depth"
oxyv.standard_name = "oxygen"

yearv.axis = "T"
mnthv.axis = "M"
lonv.axis = "X"
latv.axis = "Y"
depv.axis = "Z"

oxyv.coordinates = "t m z y x"
countv.coordinates = "t m z y x"

months = np.arange(1,13,1)
lon[lon<20] += 360.

yearv[:] = years
mnthv[:] = months
lonv[:] = lon
latv[:] = lat
depv[:] = dep
oxyv[:,:,:,:,:] = oxy_grid
countv[:,:,:,:,:] = count


data.close()
