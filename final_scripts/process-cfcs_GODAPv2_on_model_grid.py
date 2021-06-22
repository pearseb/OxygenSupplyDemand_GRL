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


#%% determine number of oxygen measurments with at least one of CFC-11, CFC-12 or SF6

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data')
glodap_all = pd.read_csv('GLODAPv2.2020_Merged_Master_File.csv', usecols=(3,4,8,9,14,15,17,26,56,60,72))
print(glodap_all.columns)

cfcs = glodap_all.fillna(0)['cfc11'] + glodap_all.fillna(0)['cfc12'] + glodap_all.fillna(0)['sf6']
cfcs[cfcs==0.0] = np.nan

cfcs = cfcs.dropna()


#%% get GLODAPv2 (2019) d13c data of DIC

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data')
glodap_cfc11 = pd.read_csv('GLODAPv2.2020_Merged_Master_File.csv', usecols=(3,4,8,9,14,15,17,26,56,57))
glodap_cfc12 = pd.read_csv('GLODAPv2.2020_Merged_Master_File.csv', usecols=(3,4,8,9,14,15,17,26,60,61))
glodap_sf6 = pd.read_csv('GLODAPv2.2020_Merged_Master_File.csv', usecols=(3,4,8,9,14,15,17,26,72,73))

print(glodap_cfc11.columns)
print(np.min(glodap_cfc11['year']), np.max(glodap_cfc11['year']))
print(np.min(glodap_cfc11['month']), np.max(glodap_cfc11['month']))
print(np.min(glodap_cfc11['longitude']), np.max(glodap_cfc11['longitude']))
print(np.min(glodap_cfc11['oxygen']), np.max(glodap_cfc11['oxygen']))
print(np.min(glodap_cfc11['cfc11']), np.max(glodap_cfc11['cfc11']))
print(np.min(glodap_cfc11['pcfc11']), np.max(glodap_cfc11['pcfc11']))

print(glodap_cfc12.columns)
print(np.min(glodap_cfc12['year']), np.max(glodap_cfc12['year']))
print(np.min(glodap_cfc12['month']), np.max(glodap_cfc12['month']))
print(np.min(glodap_cfc12['longitude']), np.max(glodap_cfc12['longitude']))
print(np.min(glodap_cfc12['oxygen']), np.max(glodap_cfc12['oxygen']))
print(np.min(glodap_cfc12['cfc12']), np.max(glodap_cfc12['cfc12']))
print(np.min(glodap_cfc12['pcfc12']), np.max(glodap_cfc12['pcfc12']))

print(glodap_sf6.columns)
print(np.min(glodap_sf6['year']), np.max(glodap_sf6['year']))
print(np.min(glodap_sf6['month']), np.max(glodap_sf6['month']))
print(np.min(glodap_sf6['longitude']), np.max(glodap_sf6['longitude']))
print(np.min(glodap_sf6['oxygen']), np.max(glodap_sf6['oxygen']))
print(np.min(glodap_sf6['sf6']), np.max(glodap_sf6['sf6']))
print(np.min(glodap_sf6['psf6']), np.max(glodap_sf6['psf6']))

### remove rows with nan values
glodap_cfc11 = glodap_cfc11.dropna()
glodap_cfc12 = glodap_cfc12.dropna()
glodap_sf6 = glodap_sf6.dropna()


### make all longitudes positive
glodap_cfc11['longitude'][glodap_cfc11['longitude']<0] += 360.0
glodap_cfc12['longitude'][glodap_cfc12['longitude']<0] += 360.0
glodap_sf6['longitude'][glodap_sf6['longitude']<0] += 360.0


print("Number of GLODAPv2 (2020) CFC-11 measurements =", len(glodap_cfc11['cfc11']))
print("Number of GLODAPv2 (2020) CFC-12 measurements =", len(glodap_cfc12['cfc12']))
print("Number of GLODAPv2 (2020) SF6 measurements =", len(glodap_sf6['sf6']))



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


#%% average data at each grid point of the model and do so for every year (CFC-11)

years = np.arange(1972, 2020, 1)

oxy_grid = np.zeros((len(years),12,31,180,360))
cfc11_grid = np.zeros((len(years),12,31,180,360))
pcfc11_grid = np.zeros((len(years),12,31,180,360))
count = np.zeros(np.shape(cfc11_grid))

# for every year
for yr,year in enumerate(years):
    # for each month 
    for mn,mnth in enumerate(np.arange(1,13,1)):
        # select only data in month
        tmp = glodap_cfc11.loc[glodap_cfc11['year'] == year]
        tmp = tmp.loc[tmp['month'] == mnth]
        if len(tmp['cfc11']) == 0:
            oxy_grid[yr,mn,:,:,:] = np.nan
            cfc11_grid[yr,mn,:,:,:] = np.nan
            pcfc11_grid[yr,mn,:,:,:] = np.nan
            continue
        
        print("YEAR", year, "MONTH", mnth, "DATA", len(tmp['cfc11']))
        tmp1 = tmp
        
        # for each grid cell
        for i in tqdm(np.arange(len(lon)), desc='Longitudes'):
            print(i,lon_bnds[i,0], lon_bnds[i,1])
            tmp2 = tmp1.loc[ (tmp1['longitude'] >= lon_bnds[i,0]) & (tmp1['longitude'] < lon_bnds[i,1]) ]
            if len(tmp2['cfc11']) > 0:
                for j in np.arange(len(lat)):
                    #print(j,lat_bnds[j,0], lat_bnds[j,1])
                    tmp3 = tmp2.loc[ (tmp2['latitude'] >= lat_bnds[j,0]) & (tmp2['latitude'] < lat_bnds[j,1]) ]
                    if len(tmp3['cfc11']) > 0:
                        for k in np.arange(len(dep)):
                            #print(k,dep_bnds[k,0], dep_bnds[k,1])
                            tmp4 = tmp3.loc[ (tmp3['depth'] >= dep_bnds[k,0]) & (tmp3['depth'] < dep_bnds[k,1]) ]
                            if len(tmp4['cfc11']) > 0:
                                print("Number of observations at year", year, "and month", mnth, "after location selection =", len(tmp4['cfc11']))
                                oxy_grid[yr,mn,k,j,i] = np.mean(tmp4['oxygen'])
                                cfc11_grid[yr,mn,k,j,i] = np.mean(tmp4['cfc11'])
                                pcfc11_grid[yr,mn,k,j,i] = np.mean(tmp4['pcfc11'])
                                count[yr,mn,k,j,i] = len(tmp4['cfc11'])
                            else:
                                oxy_grid[yr,mn,k,j,i] = np.nan
                                cfc11_grid[yr,mn,k,j,i] = np.nan
                                pcfc11_grid[yr,mn,k,j,i] = np.nan
                    else:    
                        oxy_grid[yr,mn,:,j,i] = np.nan
                        cfc11_grid[yr,mn,:,j,i] = np.nan
                        pcfc11_grid[yr,mn,:,j,i] = np.nan
            else:
                oxy_grid[yr,mn,:,:,i] = np.nan
                cfc11_grid[yr,mn,:,:,i] = np.nan
                pcfc11_grid[yr,mn,:,:,i] = np.nan
                
print(np.nanmin(cfc11_grid), np.nanmax(cfc11_grid))


for yr,year in enumerate(years):
    for mn,mnth in enumerate(np.arange(1,13,1)):
        print("Year =", year, "Month =", mnth, np.nanmin(cfc11_grid[yr,mn,:,:,:]), np.nanmax(cfc11_grid[yr,mn,:,:,:]))

print("Observations after regridding = ",np.ma.count(np.ma.masked_where(np.isnan(cfc11_grid), cfc11_grid)))


#%% save as numpy array

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\CFCs')
np.savez('GLODAPv2.2020_CFC11_1x1degrees.npz', oxy=oxy_grid, cfc11=cfc11_grid, pcfc11=pcfc11_grid)


#%% save as netcdf

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\CFCs')

os.remove('GLODAPv2.2020_CFC11_1x1degrees.nc')
data = nc.Dataset('GLODAPv2.2020_CFC11_1x1degrees.nc', 'w', format='NETCDF4_CLASSIC')

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
oxyv = data.createVariable('O2', np.float64, ('t', 'm', 'z', 'y', 'x'))
cfc11v = data.createVariable('CFC11', np.float64, ('t', 'm', 'z', 'y', 'x'))
pcfc11v = data.createVariable('pCFC11', np.float64, ('t', 'm', 'z', 'y', 'x'))
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
cfc11v.units = "pmol/kg"
pcfc11v.units = "ppt"
countv.units = "number of measurements"

yearv.standard_name = "year"
mnthv.standard_name = "month"
lonv.standard_name = "longitude"
latv.standard_name = "latitude"
depv.standard_name = "depth"
oxyv.standard_name = "Dissolved Oxygen"
cfc11v.standard_name = "CFC-11"
pcfc11v.standard_name = "pCFC-11"

yearv.axis = "T"
mnthv.axis = "M"
lonv.axis = "X"
latv.axis = "Y"
depv.axis = "Z"

oxyv.coordinates = "t m z y x"
cfc11v.coordinates = "t m z y x"
pcfc11v.coordinates = "t m z y x"
countv.coordinates = "t m z y x"

months = np.arange(1,13,1)
lon2 = lon*1
lon2[lon2<20] += 360.

yearv[:] = years
mnthv[:] = months
lonv[:] = lon2
latv[:] = lat
depv[:] = dep
oxyv[:,:,:,:,:] = oxy_grid
cfc11v[:,:,:,:,:] = cfc11_grid
pcfc11v[:,:,:,:,:] = pcfc11_grid
countv[:,:,:,:,:] = count


data.close()



#%% CFC-12

years = np.arange(1972, 2020, 1)

oxy_grid = np.zeros((len(years),12,31,180,360))
cfc12_grid = np.zeros((len(years),12,31,180,360))
pcfc12_grid = np.zeros((len(years),12,31,180,360))
count = np.zeros(np.shape(cfc12_grid))

# for every year
for yr,year in enumerate(years):
    # for each month 
    for mn,mnth in enumerate(np.arange(1,13,1)):
        # select only data in month
        tmp = glodap_cfc12.loc[glodap_cfc12['year'] == year]
        tmp = tmp.loc[tmp['month'] == mnth]
        if len(tmp['cfc12']) == 0:
            oxy_grid[yr,mn,:,:,:] = np.nan
            cfc12_grid[yr,mn,:,:,:] = np.nan
            pcfc12_grid[yr,mn,:,:,:] = np.nan
            continue
        
        print("YEAR", year, "MONTH", mnth, "DATA", len(tmp['cfc12']))
        tmp1 = tmp
        
        # for each grid cell
        for i in tqdm(np.arange(len(lon)), desc='Longitudes'):
            print(i,lon_bnds[i,0], lon_bnds[i,1])
            tmp2 = tmp1.loc[ (tmp1['longitude'] >= lon_bnds[i,0]) & (tmp1['longitude'] < lon_bnds[i,1]) ]
            if len(tmp2['cfc12']) > 0:
                for j in np.arange(len(lat)):
                    #print(j,lat_bnds[j,0], lat_bnds[j,1])
                    tmp3 = tmp2.loc[ (tmp2['latitude'] >= lat_bnds[j,0]) & (tmp2['latitude'] < lat_bnds[j,1]) ]
                    if len(tmp3['cfc12']) > 0:
                        for k in np.arange(len(dep)):
                            #print(k,dep_bnds[k,0], dep_bnds[k,1])
                            tmp4 = tmp3.loc[ (tmp3['depth'] >= dep_bnds[k,0]) & (tmp3['depth'] < dep_bnds[k,1]) ]
                            if len(tmp4['cfc12']) > 0:
                                print("Number of observations at year", year, "and month", mnth, "after location selection =", len(tmp4['cfc12']))
                                oxy_grid[yr,mn,k,j,i] = np.mean(tmp4['oxygen'])
                                cfc12_grid[yr,mn,k,j,i] = np.mean(tmp4['cfc12'])
                                pcfc12_grid[yr,mn,k,j,i] = np.mean(tmp4['pcfc12'])
                                count[yr,mn,k,j,i] = len(tmp4['cfc12'])
                            else:
                                oxy_grid[yr,mn,k,j,i] = np.nan
                                cfc12_grid[yr,mn,k,j,i] = np.nan
                                pcfc12_grid[yr,mn,k,j,i] = np.nan
                    else:    
                        oxy_grid[yr,mn,:,j,i] = np.nan
                        cfc12_grid[yr,mn,:,j,i] = np.nan
                        pcfc12_grid[yr,mn,:,j,i] = np.nan
            else:
                oxy_grid[yr,mn,:,:,i] = np.nan
                cfc12_grid[yr,mn,:,:,i] = np.nan
                pcfc12_grid[yr,mn,:,:,i] = np.nan
                
print(np.nanmin(cfc12_grid), np.nanmax(cfc12_grid))


for yr,year in enumerate(years):
    for mn,mnth in enumerate(np.arange(1,13,1)):
        print("Year =", year, "Month =", mnth, np.nanmin(cfc12_grid[yr,mn,:,:,:]), np.nanmax(cfc12_grid[yr,mn,:,:,:]))

print("Observations after regridding = ",np.ma.count(np.ma.masked_where(np.isnan(cfc12_grid), cfc12_grid)))


#%% save as numpy array

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\CFCs')
np.savez('GLODAPv2.2020_CFC12_1x1degrees.npz', oxy=oxy_grid, cfc12=cfc12_grid, pcfc12=pcfc12_grid)


#%% save as netcdf

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\CFCs')

os.remove('GLODAPv2.2020_CFC12_1x1degrees.nc')
data = nc.Dataset('GLODAPv2.2020_CFC12_1x1degrees.nc', 'w', format='NETCDF4_CLASSIC')

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
oxyv = data.createVariable('O2', np.float64, ('t', 'm', 'z', 'y', 'x'))
cfc12v = data.createVariable('CFC12', np.float64, ('t', 'm', 'z', 'y', 'x'))
pcfc12v = data.createVariable('pCFC12', np.float64, ('t', 'm', 'z', 'y', 'x'))
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
cfc12v.units = "pmol/kg"
pcfc12v.units = "ppt"
countv.units = "number of measurements"

yearv.standard_name = "year"
mnthv.standard_name = "month"
lonv.standard_name = "longitude"
latv.standard_name = "latitude"
depv.standard_name = "depth"
oxyv.standard_name = "Dissolved Oxygen"
cfc12v.standard_name = "CFC-12"
pcfc12v.standard_name = "pCFC-12"

yearv.axis = "T"
mnthv.axis = "M"
lonv.axis = "X"
latv.axis = "Y"
depv.axis = "Z"

oxyv.coordinates = "t m z y x"
cfc12v.coordinates = "t m z y x"
pcfc12v.coordinates = "t m z y x"
countv.coordinates = "t m z y x"

months = np.arange(1,13,1)
lon2 = lon*1
lon2[lon2<20] += 360.

yearv[:] = years
mnthv[:] = months
lonv[:] = lon2
latv[:] = lat
depv[:] = dep
oxyv[:,:,:,:,:] = oxy_grid
cfc12v[:,:,:,:,:] = cfc12_grid
pcfc12v[:,:,:,:,:] = pcfc12_grid
countv[:,:,:,:,:] = count


data.close()


#%% SF6

years = np.arange(1972, 2020, 1)

oxy_grid = np.zeros((len(years),12,31,180,360))
sf6_grid = np.zeros((len(years),12,31,180,360))
psf6_grid = np.zeros((len(years),12,31,180,360))
count = np.zeros(np.shape(sf6_grid))

# for every year
for yr,year in enumerate(years):
    # for each month 
    for mn,mnth in enumerate(np.arange(1,13,1)):
        # select only data in month
        tmp = glodap_sf6.loc[glodap_sf6['year'] == year]
        tmp = tmp.loc[tmp['month'] == mnth]
        if len(tmp['sf6']) == 0:
            oxy_grid[yr,mn,:,:,:] = np.nan
            sf6_grid[yr,mn,:,:,:] = np.nan
            psf6_grid[yr,mn,:,:,:] = np.nan
            continue
        
        print("YEAR", year, "MONTH", mnth, "DATA", len(tmp['sf6']))
        tmp1 = tmp
        
        # for each grid cell
        for i in tqdm(np.arange(len(lon)), desc='Longitudes'):
            print(i,lon_bnds[i,0], lon_bnds[i,1])
            tmp2 = tmp1.loc[ (tmp1['longitude'] >= lon_bnds[i,0]) & (tmp1['longitude'] < lon_bnds[i,1]) ]
            if len(tmp2['sf6']) > 0:
                for j in np.arange(len(lat)):
                    #print(j,lat_bnds[j,0], lat_bnds[j,1])
                    tmp3 = tmp2.loc[ (tmp2['latitude'] >= lat_bnds[j,0]) & (tmp2['latitude'] < lat_bnds[j,1]) ]
                    if len(tmp3['sf6']) > 0:
                        for k in np.arange(len(dep)):
                            #print(k,dep_bnds[k,0], dep_bnds[k,1])
                            tmp4 = tmp3.loc[ (tmp3['depth'] >= dep_bnds[k,0]) & (tmp3['depth'] < dep_bnds[k,1]) ]
                            if len(tmp4['sf6']) > 0:
                                print("Number of observations at year", year, "and month", mnth, "after location selection =", len(tmp4['sf6']))
                                oxy_grid[yr,mn,k,j,i] = np.mean(tmp4['oxygen'])
                                sf6_grid[yr,mn,k,j,i] = np.mean(tmp4['sf6'])
                                psf6_grid[yr,mn,k,j,i] = np.mean(tmp4['psf6'])
                                count[yr,mn,k,j,i] = len(tmp4['sf6'])
                            else:
                                oxy_grid[yr,mn,k,j,i] = np.nan
                                sf6_grid[yr,mn,k,j,i] = np.nan
                                psf6_grid[yr,mn,k,j,i] = np.nan
                    else:    
                        oxy_grid[yr,mn,:,j,i] = np.nan
                        sf6_grid[yr,mn,:,j,i] = np.nan
                        psf6_grid[yr,mn,:,j,i] = np.nan
            else:
                oxy_grid[yr,mn,:,:,i] = np.nan
                sf6_grid[yr,mn,:,:,i] = np.nan
                psf6_grid[yr,mn,:,:,i] = np.nan
                
print(np.nanmin(sf6_grid), np.nanmax(sf6_grid))


for yr,year in enumerate(years):
    for mn,mnth in enumerate(np.arange(1,13,1)):
        print("Year =", year, "Month =", mnth, np.nanmin(sf6_grid[yr,mn,:,:,:]), np.nanmax(sf6_grid[yr,mn,:,:,:]))

print("Observations after regridding = ",np.ma.count(np.ma.masked_where(np.isnan(sf6_grid), sf6_grid)))


#%% save as numpy array

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\CFCs')
np.savez('GLODAPv2.2020_SF6_1x1degrees.npz', oxy=oxy_grid, sf6=sf6_grid, psf6=psf6_grid)


#%% save as netcdf

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\Data\\CFCs')

os.remove('GLODAPv2.2020_SF6_1x1degrees.nc')
data = nc.Dataset('GLODAPv2.2020_SF6_1x1degrees.nc', 'w', format='NETCDF4_CLASSIC')

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
oxyv = data.createVariable('O2', np.float64, ('t', 'm', 'z', 'y', 'x'))
sf6v = data.createVariable('SF6', np.float64, ('t', 'm', 'z', 'y', 'x'))
psf6v = data.createVariable('pSF6', np.float64, ('t', 'm', 'z', 'y', 'x'))
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
sf6v.units = "fmol/kg"
psf6v.units = "ppt"
countv.units = "number of measurements"

yearv.standard_name = "year"
mnthv.standard_name = "month"
lonv.standard_name = "longitude"
latv.standard_name = "latitude"
depv.standard_name = "depth"
oxyv.standard_name = "Dissolved Oxygen"
sf6v.standard_name = "SF6"
psf6v.standard_name = "pSF6"

yearv.axis = "T"
mnthv.axis = "M"
lonv.axis = "X"
latv.axis = "Y"
depv.axis = "Z"

oxyv.coordinates = "t m z y x"
sf6v.coordinates = "t m z y x"
psf6v.coordinates = "t m z y x"
countv.coordinates = "t m z y x"

months = np.arange(1,13,1)
lon2 = lon*1
lon2[lon2<20] += 360.

yearv[:] = years
mnthv[:] = months
lonv[:] = lon2
latv[:] = lat
depv[:] = dep
oxyv[:,:,:,:,:] = oxy_grid
sf6v[:,:,:,:,:] = sf6_grid
psf6v[:,:,:,:,:] = psf6_grid
countv[:,:,:,:,:] = count


data.close()
