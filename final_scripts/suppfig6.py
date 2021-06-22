# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:45:34 2021

Purpose
-------
    Figure of linear oxygen, TOU and ideal age trends averaged over 26.5-27.0 and
    27.0-27.4 potential density bounds 
    (hatching represents where the signal exceeds interannual noise)
    
    
@author: pearseb
"""

#%% imports

from __future__ import unicode_literals

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cmocean
import cmocean.cm as cmo
import seaborn as sb
sb.set(style='ticks')

from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy
import cartopy.crs as ccrs

# for colouring scatter plot according to density of data
from scipy.stats import gaussian_kde

# for lowess smoothing
from statsmodels.nonparametric.smoothers_lowess import lowess


#%% get data

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")


data = nc.Dataset("suppfig_ipsl_age.nc", 'r')
aou_ipsl_tre = np.squeeze(data.variables['AOU_TRE_ZAVE'][...])
wma_ipsl_tre = np.squeeze(data.variables['WMA_TRE_ZAVE'][...])

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()


p18 = np.genfromtxt('stations_P18.txt', delimiter='\t')
p16n = np.genfromtxt('stations_P16N.txt', delimiter='\t')
a135 = np.genfromtxt('stations_A13.5.txt', delimiter='\t')
a16s = np.genfromtxt('stations_A16S.txt', delimiter='\t')
a16s[:,0] = a16s[:,0] + a16s[:,1]/60.0 
a16s[:,1] = a16s[:,2] + a16s[:,3]/60.0 
a16s = a16s[:,0:2]
a16n = np.genfromtxt('stations_A16N.txt', delimiter='\t')
a16n[:,0] = a16n[:,0] + a16n[:,1]/60.0 
a16n[:,1] = a16n[:,2] + a16n[:,3]/60.0 
a16n = a16n[:,0:2]
a10 = np.genfromtxt('stations_A10.txt', delimiter='\t')
a10[:,0] = a10[:,0] + a10[:,1]/60.0 
a10[:,1] = a10[:,2] + a10[:,3]/60.0 
a10 = a10[:,0:2]
i05 = np.genfromtxt('stations_I05.txt', delimiter='\t')
i05[:,0] = i05[:,0] + i05[:,1]/60.0 
i05[:,1] = i05[:,2] + i05[:,3]/60.0 
i05 = i05[:,0:2]



#%% get data for A16S

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_A16S.nc', 'r')
oxy_m_A16S = np.squeeze(data.variables['AOU_M_DIF'][...])/25*10.
oxy_o_A16S = np.squeeze(data.variables['AOU_O_DIF'][...])/25*10.
wma_m_A16S = np.squeeze(data.variables['WMA_M_DIF'][...])/25*10.
wma_o_A16S = np.squeeze(data.variables['WMA_O_DIF'][...])/25*10.
oxy_m_A16S_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/25*10.
oxy_o_A16S_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/25*10.
wma_m_A16S_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/25*10.
wma_o_A16S_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/25*10.
rho_m_2014_A16S = np.squeeze(data.variables['RHO_M_2014'][...])
rho_o_2014_A16S = np.squeeze(data.variables['RHO_O_2014'][...])
rho_m_1989_A16S = np.squeeze(data.variables['RHO_M_1989'][...])
rho_o_1989_A16S = np.squeeze(data.variables['RHO_O_1989'][...])
oxy_m_rho_A16S = np.squeeze(data.variables['OXY_M_DIF_26P5_27P4'][...])/25*10.
oxy_o_rho_A16S = np.squeeze(data.variables['OXY_O_DIF_26P5_27P4'][...])/25*10.
sat_m_rho_A16S = np.squeeze(data.variables['SAT_M_DIF_26P5_27P4'][...])/25*10.
sat_o_rho_A16S = np.squeeze(data.variables['SAT_O_DIF_26P5_27P4'][...])/25*10.
aou_m_rho_A16S = np.squeeze(data.variables['AOU_M_DIF_26P5_27P4'][...])/25*10.*(-1)
aou_o_rho_A16S = np.squeeze(data.variables['AOU_O_DIF_26P5_27P4'][...])/25*10.*(-1)
age_m_rho_A16S = np.squeeze(data.variables['WMA_M_DIF_26P5_27P4'][...])/25*10.
age_o_rho_A16S = np.squeeze(data.variables['WMA_O_DIF_26P5_27P4'][...])/25*10.
dem_m_rho_A16S = np.squeeze(data.variables['DEM_M_DIF_26P5_27P4'][...])/25*10.
dem_o_rho_A16S = np.squeeze(data.variables['DEM_O_DIF_26P5_27P4'][...])/25*10.

sup_m_rho_A16S = aou_m_rho_A16S - dem_m_rho_A16S
sup_o_rho_A16S = aou_o_rho_A16S - dem_o_rho_A16S
fdem_m_rho_A16S = dem_m_rho_A16S/aou_m_rho_A16S
fsup_m_rho_A16S = sup_m_rho_A16S/aou_m_rho_A16S
fdem_o_rho_A16S = dem_o_rho_A16S/aou_o_rho_A16S
fsup_o_rho_A16S = sup_o_rho_A16S/aou_o_rho_A16S


dep_A16S = data.variables['DEPTHT11_23'][...]
lat_A16S = data.variables['ETOPO60Y121_100'][...]
lon_A16S = np.zeros(len(lat_A16S))
for i,val in enumerate(lat_A16S):
    tmp = np.abs(val - a16s[:,0])
    lon_A16S[i] = a16s[int(np.where(tmp==np.min(tmp))[0][0]),1]

oxy_m_A16S_zave_nans = np.ma.getdata(oxy_m_A16S_zave)
oxy_m_A16S_zave_nans[oxy_m_A16S_zave_nans<-10000] = np.nan
wma_m_A16S_zave_nans = np.ma.getdata(wma_m_A16S_zave)
wma_m_A16S_zave_nans[wma_m_A16S_zave_nans<-10000] = np.nan
oxy_m_A16S_zave_smooth = lowess(oxy_m_A16S_zave_nans, lat_A16S, frac=0.2, delta=5)
wma_m_A16S_zave_smooth = lowess(wma_m_A16S_zave_nans, lat_A16S, frac=0.2, delta=5)

oxy_m_A16S_section = np.zeros((len(oxy_m_A16S_zave_smooth[:,0]),len(oxy_m_A16S_zave_smooth[0,:])+1)) 
oxy_m_A16S_section[:,1] = oxy_m_A16S_zave_smooth[:,0]
oxy_m_A16S_section[:,2] = oxy_m_A16S_zave_smooth[:,1]
for i,val in enumerate(oxy_m_A16S_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16S)
    oxy_m_A16S_section[i,0] = lon_A16S[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_A16S_section = np.zeros((len(wma_m_A16S_zave_smooth[:,0]),len(wma_m_A16S_zave_smooth[0,:])+1)) 
wma_m_A16S_section[:,1] = wma_m_A16S_zave_smooth[:,0]
wma_m_A16S_section[:,2] = wma_m_A16S_zave_smooth[:,1]
for i,val in enumerate(wma_m_A16S_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16S)
    wma_m_A16S_section[i,0] = lon_A16S[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_A16S_zave_nans = np.ma.getdata(oxy_o_A16S_zave)
oxy_o_A16S_zave_nans[oxy_o_A16S_zave_nans<-10000] = np.nan
wma_o_A16S_zave_nans = np.ma.getdata(wma_o_A16S_zave)
wma_o_A16S_zave_nans[wma_o_A16S_zave_nans<-10000] = np.nan
oxy_o_A16S_zave_smooth = lowess(oxy_o_A16S_zave_nans, lat_A16S, frac=0.2, delta=5)
wma_o_A16S_zave_smooth = lowess(wma_o_A16S_zave_nans, lat_A16S, frac=0.2, delta=5)

oxy_o_A16S_section = np.zeros((len(oxy_o_A16S_zave_smooth[:,0]),len(oxy_o_A16S_zave_smooth[0,:])+1)) 
oxy_o_A16S_section[:,1] = oxy_o_A16S_zave_smooth[:,0]
oxy_o_A16S_section[:,2] = oxy_o_A16S_zave_smooth[:,1]
for i,val in enumerate(oxy_o_A16S_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16S)
    oxy_o_A16S_section[i,0] = lon_A16S[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_A16S_section = np.zeros((len(wma_o_A16S_zave_smooth[:,0]),len(wma_o_A16S_zave_smooth[0,:])+1)) 
wma_o_A16S_section[:,1] = wma_o_A16S_zave_smooth[:,0]
wma_o_A16S_section[:,2] = wma_o_A16S_zave_smooth[:,1]
for i,val in enumerate(wma_o_A16S_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16S)
    wma_o_A16S_section[i,0] = lon_A16S[int(np.where(tmp==np.min(tmp))[0][0])]


data.close()


#%% get data for A16N

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_A16N.nc', 'r')
oxy_m_A16N = np.squeeze(data.variables['AOU_M_DIF'][...])/25*10.
oxy_o_A16N = np.squeeze(data.variables['AOU_O_DIF'][...])/25*10.
wma_m_A16N = np.squeeze(data.variables['WMA_M_DIF'][...])/25*10.
wma_o_A16N = np.squeeze(data.variables['WMA_O_DIF'][...])/25*10.
oxy_m_A16N_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/25*10.
oxy_o_A16N_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/25*10.
wma_m_A16N_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/25*10.
wma_o_A16N_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/25*10.
rho_m_2013_A16N = np.squeeze(data.variables['RHO_M_2013'][...])
rho_o_2013_A16N = np.squeeze(data.variables['RHO_O_2013'][...])
rho_m_1988_A16N = np.squeeze(data.variables['RHO_M_1988'][...])
rho_o_1988_A16N = np.squeeze(data.variables['RHO_O_1988'][...])
oxy_m_rho_A16N = np.squeeze(data.variables['OXY_M_DIF_27P0_27P6'][...])/25*10.
oxy_o_rho_A16N = np.squeeze(data.variables['OXY_O_DIF_27P0_27P6'][...])/25*10.
sat_m_rho_A16N = np.squeeze(data.variables['SAT_M_DIF_27P0_27P6'][...])/25*10.
sat_o_rho_A16N = np.squeeze(data.variables['SAT_O_DIF_27P0_27P6'][...])/25*10.
aou_m_rho_A16N = np.squeeze(data.variables['AOU_M_DIF_27P0_27P6'][...])/25*10.*(-1)
aou_o_rho_A16N = np.squeeze(data.variables['AOU_O_DIF_27P0_27P6'][...])/25*10.*(-1)
age_m_rho_A16N = np.squeeze(data.variables['WMA_M_DIF_27P0_27P6'][...])/25*10.
age_o_rho_A16N = np.squeeze(data.variables['WMA_O_DIF_27P0_27P6'][...])/25*10.
dem_m_rho_A16N = np.squeeze(data.variables['DEM_M_DIF_27P0_27P6'][...])/25*10.
dem_o_rho_A16N = np.squeeze(data.variables['DEM_O_DIF_27P0_27P6'][...])/25*10.

sup_m_rho_A16N = aou_m_rho_A16N - dem_m_rho_A16N
sup_o_rho_A16N = aou_o_rho_A16N - dem_o_rho_A16N
fdem_m_rho_A16N = dem_m_rho_A16N/aou_m_rho_A16N
fsup_m_rho_A16N = sup_m_rho_A16N/aou_m_rho_A16N
fdem_o_rho_A16N = dem_o_rho_A16N/aou_o_rho_A16N
fsup_o_rho_A16N = sup_o_rho_A16N/aou_o_rho_A16N


dep_A16N = data.variables['DEPTHT11_23'][...]
lat_A16N = data.variables['ETOPO60Y181_160'][...]
lon_A16N = np.zeros(len(lat_A16N))
for i,val in enumerate(lat_A16N):
    tmp = np.abs(val - a16n[:,0])
    lon_A16N[i] = a16n[int(np.where(tmp==np.min(tmp))[0][0]),1]

# smoothing
oxy_m_A16N_zave_nans = np.ma.getdata(oxy_m_A16N_zave)
oxy_m_A16N_zave_nans[oxy_m_A16N_zave_nans<-10000] = np.nan
wma_m_A16N_zave_nans = np.ma.getdata(wma_m_A16N_zave)
wma_m_A16N_zave_nans[wma_m_A16N_zave_nans<-10000] = np.nan
oxy_m_A16N_zave_smooth = lowess(oxy_m_A16N_zave_nans, lat_A16N, frac=0.2, delta=5)
wma_m_A16N_zave_smooth = lowess(wma_m_A16N_zave_nans, lat_A16N, frac=0.2, delta=5)

oxy_m_A16N_section = np.zeros((len(oxy_m_A16N_zave_smooth[:,0]),len(oxy_m_A16N_zave_smooth[0,:])+1)) 
oxy_m_A16N_section[:,1] = oxy_m_A16N_zave_smooth[:,0]
oxy_m_A16N_section[:,2] = oxy_m_A16N_zave_smooth[:,1]
for i,val in enumerate(oxy_m_A16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16N)
    oxy_m_A16N_section[i,0] = lon_A16N[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_A16N_section = np.zeros((len(wma_m_A16N_zave_smooth[:,0]),len(wma_m_A16N_zave_smooth[0,:])+1)) 
wma_m_A16N_section[:,1] = wma_m_A16N_zave_smooth[:,0]
wma_m_A16N_section[:,2] = wma_m_A16N_zave_smooth[:,1]
for i,val in enumerate(wma_m_A16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16N)
    wma_m_A16N_section[i,0] = lon_A16N[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_A16N_zave_nans = np.ma.getdata(oxy_o_A16N_zave)
oxy_o_A16N_zave_nans[oxy_o_A16N_zave_nans<-10000] = np.nan
wma_o_A16N_zave_nans = np.ma.getdata(wma_o_A16N_zave)
wma_o_A16N_zave_nans[wma_o_A16N_zave_nans<-10000] = np.nan
oxy_o_A16N_zave_smooth = lowess(oxy_o_A16N_zave_nans, lat_A16N, frac=0.2, delta=5)
wma_o_A16N_zave_smooth = lowess(wma_o_A16N_zave_nans, lat_A16N, frac=0.2, delta=5)

oxy_o_A16N_section = np.zeros((len(oxy_o_A16N_zave_smooth[:,0]),len(oxy_o_A16N_zave_smooth[0,:])+1)) 
oxy_o_A16N_section[:,1] = oxy_o_A16N_zave_smooth[:,0]
oxy_o_A16N_section[:,2] = oxy_o_A16N_zave_smooth[:,1]
for i,val in enumerate(oxy_o_A16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16N)
    oxy_o_A16N_section[i,0] = lon_A16N[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_A16N_section = np.zeros((len(wma_o_A16N_zave_smooth[:,0]),len(wma_o_A16N_zave_smooth[0,:])+1)) 
wma_o_A16N_section[:,1] = wma_o_A16N_zave_smooth[:,0]
wma_o_A16N_section[:,2] = wma_o_A16N_zave_smooth[:,1]
for i,val in enumerate(wma_o_A16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A16N)
    wma_o_A16N_section[i,0] = lon_A16N[int(np.where(tmp==np.min(tmp))[0][0])]

data.close()


#%% get data for P16N


os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_P16N.nc', 'r')
oxy_m_P16N = np.squeeze(data.variables['AOU_M_DIF'][...])/7*10.
oxy_o_P16N = np.squeeze(data.variables['AOU_O_DIF'][...])/7*10.
wma_m_P16N = np.squeeze(data.variables['WMA_M_DIF'][...])/7*10.
wma_o_P16N = np.squeeze(data.variables['WMA_O_DIF'][...])/7*10.
oxy_m_P16N_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/7*10.
oxy_o_P16N_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/7*10.
wma_m_P16N_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/7*10.
wma_o_P16N_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/7*10.
rho_m_2015_P16N = np.squeeze(data.variables['RHO_M_2015'][...])
rho_o_2015_P16N = np.squeeze(data.variables['RHO_O_2015'][...])
rho_m_2008_P16N = np.squeeze(data.variables['RHO_M_2008'][...])
rho_o_2008_P16N = np.squeeze(data.variables['RHO_O_2008'][...])
oxy_m_rho_P16N = np.squeeze(data.variables['OXY_M_DIF_26P0_26P5'][...])/7*10.
oxy_o_rho_P16N = np.squeeze(data.variables['OXY_O_DIF_26P0_26P5'][...])/7*10.
sat_m_rho_P16N = np.squeeze(data.variables['SAT_M_DIF_26P0_26P5'][...])/7*10.
sat_o_rho_P16N = np.squeeze(data.variables['SAT_O_DIF_26P0_26P5'][...])/7*10.
aou_m_rho_P16N = np.squeeze(data.variables['AOU_M_DIF_26P0_26P5'][...])/7*10.*(-1)
aou_o_rho_P16N = np.squeeze(data.variables['AOU_O_DIF_26P0_26P5'][...])/7*10.*(-1)
age_m_rho_P16N = np.squeeze(data.variables['WMA_M_DIF_26P0_26P5'][...])/7*10.
age_o_rho_P16N = np.squeeze(data.variables['WMA_O_DIF_26P0_26P5'][...])/7*10.
dem_m_rho_P16N = np.squeeze(data.variables['DEM_M_DIF_26P0_26P5'][...])/7*10.
dem_o_rho_P16N = np.squeeze(data.variables['DEM_O_DIF_26P0_26P5'][...])/7*10.

sup_m_rho_P16N = aou_m_rho_P16N - dem_m_rho_P16N
sup_o_rho_P16N = aou_o_rho_P16N - dem_o_rho_P16N
fdem_m_rho_P16N = dem_m_rho_P16N/aou_m_rho_P16N
fsup_m_rho_P16N = sup_m_rho_P16N/aou_m_rho_P16N
fdem_o_rho_P16N = dem_o_rho_P16N/aou_o_rho_P16N
fsup_o_rho_P16N = sup_o_rho_P16N/aou_o_rho_P16N


dep_P16N = data.variables['DEPTHT11_22'][...]
lat_P16N = data.variables['ETOPO60Y171_150'][...]
lon_P16N = np.zeros(len(lat_P16N))
for i,val in enumerate(lat_P16N):
    tmp = np.abs(val - p16n[:,0])
    lon_P16N[i] = p16n[int(np.where(tmp==np.min(tmp))[0][0]),1]

# smoothing
oxy_m_P16N_zave_nans = np.ma.getdata(oxy_m_P16N_zave)
oxy_m_P16N_zave_nans[oxy_m_P16N_zave_nans<-10000] = np.nan
wma_m_P16N_zave_nans = np.ma.getdata(wma_m_P16N_zave)
wma_m_P16N_zave_nans[wma_m_P16N_zave_nans<-10000] = np.nan
oxy_m_P16N_zave_smooth = lowess(oxy_m_P16N_zave_nans, lat_P16N, frac=0.2, delta=5)
wma_m_P16N_zave_smooth = lowess(wma_m_P16N_zave_nans, lat_P16N, frac=0.2, delta=5)

oxy_m_P16N_section = np.zeros((len(oxy_m_P16N_zave_smooth[:,0]),len(oxy_m_P16N_zave_smooth[0,:])+1)) 
oxy_m_P16N_section[:,1] = oxy_m_P16N_zave_smooth[:,0]
oxy_m_P16N_section[:,2] = oxy_m_P16N_zave_smooth[:,1]
for i,val in enumerate(oxy_m_P16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P16N)
    oxy_m_P16N_section[i,0] = lon_P16N[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_P16N_section = np.zeros((len(wma_m_P16N_zave_smooth[:,0]),len(wma_m_P16N_zave_smooth[0,:])+1)) 
wma_m_P16N_section[:,1] = wma_m_P16N_zave_smooth[:,0]
wma_m_P16N_section[:,2] = wma_m_P16N_zave_smooth[:,1]
for i,val in enumerate(wma_m_P16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P16N)
    wma_m_P16N_section[i,0] = lon_P16N[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_P16N_zave_nans = np.ma.getdata(oxy_o_P16N_zave)
oxy_o_P16N_zave_nans[oxy_o_P16N_zave_nans<-10000] = np.nan
wma_o_P16N_zave_nans = np.ma.getdata(wma_o_P16N_zave)
wma_o_P16N_zave_nans[wma_o_P16N_zave_nans<-10000] = np.nan
oxy_o_P16N_zave_smooth = lowess(oxy_o_P16N_zave_nans, lat_P16N, frac=0.2, delta=5)
wma_o_P16N_zave_smooth = lowess(wma_o_P16N_zave_nans, lat_P16N, frac=0.2, delta=5)

oxy_o_P16N_section = np.zeros((len(oxy_o_P16N_zave_smooth[:,0]),len(oxy_o_P16N_zave_smooth[0,:])+1)) 
oxy_o_P16N_section[:,1] = oxy_o_P16N_zave_smooth[:,0]
oxy_o_P16N_section[:,2] = oxy_o_P16N_zave_smooth[:,1]
for i,val in enumerate(oxy_o_P16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P16N)
    oxy_o_P16N_section[i,0] = lon_P16N[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_P16N_section = np.zeros((len(wma_o_P16N_zave_smooth[:,0]),len(wma_o_P16N_zave_smooth[0,:])+1)) 
wma_o_P16N_section[:,1] = wma_o_P16N_zave_smooth[:,0]
wma_o_P16N_section[:,2] = wma_o_P16N_zave_smooth[:,1]
for i,val in enumerate(wma_o_P16N_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P16N)
    wma_o_P16N_section[i,0] = lon_P16N[int(np.where(tmp==np.min(tmp))[0][0])]

data.close()


#%% get data for P18

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_P18.nc', 'r')
oxy_m_P18 = np.squeeze(data.variables['AOU_M_DIF'][...])/9*10.
oxy_o_P18 = np.squeeze(data.variables['AOU_O_DIF'][...])/9*10.
wma_m_P18 = np.squeeze(data.variables['WMA_M_DIF'][...])/9*10.
wma_o_P18 = np.squeeze(data.variables['WMA_O_DIF'][...])/9*10.
oxy_m_P18_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/9*10.
oxy_o_P18_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/9*10.
wma_m_P18_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/9*10.
wma_o_P18_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/9*10.
rho_m_2016_2019_P18 = np.squeeze(data.variables['RHO_M_2016_2017'][...])
rho_o_2016_2019_P18 = np.squeeze(data.variables['RHO_O_2016_2017'][...])
rho_m_2009_2008_P18 = np.squeeze(data.variables['RHO_M_2007_2008'][...])
rho_o_2009_2008_P18 = np.squeeze(data.variables['RHO_O_2007_2008'][...])
oxy_m_rho_P18 = np.squeeze(data.variables['OXY_M_DIF_26P5_27P0'][...])/9*10.
oxy_o_rho_P18 = np.squeeze(data.variables['OXY_O_DIF_26P5_27P0'][...])/9*10.
sat_m_rho_P18 = np.squeeze(data.variables['SAT_M_DIF_26P5_27P0'][...])/9*10.
sat_o_rho_P18 = np.squeeze(data.variables['SAT_O_DIF_26P5_27P0'][...])/9*10.
aou_m_rho_P18 = np.squeeze(data.variables['AOU_M_DIF_26P5_27P0'][...])/9*10.*(-1)
aou_o_rho_P18 = np.squeeze(data.variables['AOU_O_DIF_26P5_27P0'][...])/9*10.*(-1)
age_m_rho_P18 = np.squeeze(data.variables['WMA_M_DIF_26P5_27P0'][...])/9*10.
age_o_rho_P18 = np.squeeze(data.variables['WMA_O_DIF_26P5_27P0'][...])/9*10.
dem_m_rho_P18 = np.squeeze(data.variables['DEM_M_DIF_26P5_27P0'][...])/9*10.
dem_o_rho_P18 = np.squeeze(data.variables['DEM_O_DIF_26P5_27P0'][...])/9*10.

sup_m_rho_P18 = aou_m_rho_P18 - dem_m_rho_P18
sup_o_rho_P18 = aou_o_rho_P18 - dem_o_rho_P18
fdem_m_rho_P18 = dem_m_rho_P18/aou_m_rho_P18
fsup_m_rho_P18 = sup_m_rho_P18/aou_m_rho_P18
fdem_o_rho_P18 = dem_o_rho_P18/aou_o_rho_P18
fsup_o_rho_P18 = sup_o_rho_P18/aou_o_rho_P18


dep_P18 = data.variables['DEPTHT11_22'][...]
lat_P18 = data.variables['ETOPO60Y111_70'][...]
lon_P18 = np.zeros(len(lat_P18))
for i,val in enumerate(lat_P18):
    tmp = np.abs(val - p18[:,0])
    lon_P18[i] = p18[int(np.where(tmp==np.min(tmp))[0][0]),1]

# smoothing
oxy_m_P18_zave_nans = np.ma.getdata(oxy_m_P18_zave)
oxy_m_P18_zave_nans[oxy_m_P18_zave_nans<-10000] = np.nan
wma_m_P18_zave_nans = np.ma.getdata(wma_m_P18_zave)
wma_m_P18_zave_nans[wma_m_P18_zave_nans<-10000] = np.nan
oxy_m_P18_zave_smooth = lowess(oxy_m_P18_zave_nans, lat_P18, frac=0.2, delta=5)
wma_m_P18_zave_smooth = lowess(wma_m_P18_zave_nans, lat_P18, frac=0.2, delta=5)

oxy_m_P18_section = np.zeros((len(oxy_m_P18_zave_smooth[:,0]),len(oxy_m_P18_zave_smooth[0,:])+1)) 
oxy_m_P18_section[:,1] = oxy_m_P18_zave_smooth[:,0]
oxy_m_P18_section[:,2] = oxy_m_P18_zave_smooth[:,1]
for i,val in enumerate(oxy_m_P18_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P18)
    oxy_m_P18_section[i,0] = lon_P18[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_P18_section = np.zeros((len(wma_m_P18_zave_smooth[:,0]),len(wma_m_P18_zave_smooth[0,:])+1)) 
wma_m_P18_section[:,1] = wma_m_P18_zave_smooth[:,0]
wma_m_P18_section[:,2] = wma_m_P18_zave_smooth[:,1]
for i,val in enumerate(wma_m_P18_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P18)
    wma_m_P18_section[i,0] = lon_P18[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_P18_zave_nans = np.ma.getdata(oxy_o_P18_zave)
oxy_o_P18_zave_nans[oxy_o_P18_zave_nans<-10000] = np.nan
wma_o_P18_zave_nans = np.ma.getdata(wma_o_P18_zave)
wma_o_P18_zave_nans[wma_o_P18_zave_nans<-10000] = np.nan
oxy_o_P18_zave_smooth = lowess(oxy_o_P18_zave_nans, lat_P18, frac=0.2, delta=5)
wma_o_P18_zave_smooth = lowess(wma_o_P18_zave_nans, lat_P18, frac=0.2, delta=5)

oxy_o_P18_section = np.zeros((len(oxy_o_P18_zave_smooth[:,0]),len(oxy_o_P18_zave_smooth[0,:])+1)) 
oxy_o_P18_section[:,1] = oxy_o_P18_zave_smooth[:,0]
oxy_o_P18_section[:,2] = oxy_o_P18_zave_smooth[:,1]
for i,val in enumerate(oxy_o_P18_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P18)
    oxy_o_P18_section[i,0] = lon_P18[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_P18_section = np.zeros((len(wma_o_P18_zave_smooth[:,0]),len(wma_o_P18_zave_smooth[0,:])+1)) 
wma_o_P18_section[:,1] = wma_o_P18_zave_smooth[:,0]
wma_o_P18_section[:,2] = wma_o_P18_zave_smooth[:,1]
for i,val in enumerate(wma_o_P18_zave_smooth[:,0]):
    tmp = np.abs(val - lat_P18)
    wma_o_P18_section[i,0] = lon_P18[int(np.where(tmp==np.min(tmp))[0][0])]

data.close()




#%% get data for I05

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_I05.nc', 'r')
oxy_m_I05 = np.squeeze(data.variables['AOU_M_DIF'][...])/22*10.
oxy_o_I05 = np.squeeze(data.variables['AOU_O_DIF'][...])/22*10.
wma_m_I05 = np.squeeze(data.variables['WMA_M_DIF'][...])/22*10.
wma_o_I05 = np.squeeze(data.variables['WMA_O_DIF'][...])/22*10.
oxy_m_I05_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/22*10.
oxy_o_I05_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/22*10.
wma_m_I05_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/22*10.
wma_o_I05_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/22*10.
rho_m_2009_I05 = np.squeeze(data.variables['RHO_M_2009'][...])
rho_o_2009_I05 = np.squeeze(data.variables['RHO_O_2009'][...])
rho_m_1987_I05 = np.squeeze(data.variables['RHO_M_1987'][...])
rho_o_1987_I05 = np.squeeze(data.variables['RHO_O_1987'][...])
oxy_m_rho_I05 = np.squeeze(data.variables['OXY_M_DIF_26P5_27P0'][...])/22*10.
oxy_o_rho_I05 = np.squeeze(data.variables['OXY_O_DIF_26P5_27P0'][...])/22*10.
sat_m_rho_I05 = np.squeeze(data.variables['SAT_M_DIF_26P5_27P0'][...])/22*10.
sat_o_rho_I05 = np.squeeze(data.variables['SAT_O_DIF_26P5_27P0'][...])/22*10.
aou_m_rho_I05 = np.squeeze(data.variables['AOU_M_DIF_26P5_27P0'][...])/22*10.*(-1)
aou_o_rho_I05 = np.squeeze(data.variables['AOU_O_DIF_26P5_27P0'][...])/22*10.*(-1)
age_m_rho_I05 = np.squeeze(data.variables['WMA_M_DIF_26P5_27P0'][...])/22*10.
age_o_rho_I05 = np.squeeze(data.variables['WMA_O_DIF_26P5_27P0'][...])/22*10.
dem_m_rho_I05 = np.squeeze(data.variables['DEM_M_DIF_26P5_27P0'][...])/22*10.
dem_o_rho_I05 = np.squeeze(data.variables['DEM_O_DIF_26P5_27P0'][...])/22*10.

sup_m_rho_I05 = aou_m_rho_I05 - dem_m_rho_I05
sup_o_rho_I05 = aou_o_rho_I05 - dem_o_rho_I05
fdem_m_rho_I05 = dem_m_rho_I05/aou_m_rho_I05
fsup_m_rho_I05 = sup_m_rho_I05/aou_m_rho_I05
fdem_o_rho_I05 = dem_o_rho_I05/aou_o_rho_I05
fsup_o_rho_I05 = sup_o_rho_I05/aou_o_rho_I05


dep_I05 = data.variables['DEPTHT11_22'][...]
lon_I05 = data.variables['ETOPO60X16_100'][...]
lat_I05 = np.zeros(len(lon_I05))
for i,val in enumerate(lon_I05):
    tmp = np.abs(val - i05[:,1])
    lat_I05[i] = i05[int(np.where(tmp==np.min(tmp))[0][0]),0]


# smoothing
oxy_m_I05_zave_nans = np.ma.getdata(oxy_m_I05_zave)
oxy_m_I05_zave_nans[oxy_m_I05_zave_nans<-10000] = np.nan
oxy_m_I05_zave_nans[oxy_m_I05_zave_nans<-50] = np.nan
wma_m_I05_zave_nans = np.ma.getdata(wma_m_I05_zave)
wma_m_I05_zave_nans[wma_m_I05_zave_nans<-10000] = np.nan
oxy_m_I05_zave_smooth = lowess(oxy_m_I05_zave_nans, lon_I05, frac=0.2, delta=5)
wma_m_I05_zave_smooth = lowess(wma_m_I05_zave_nans, lon_I05, frac=0.2, delta=5)

oxy_m_I05_section = np.zeros((len(oxy_m_I05_zave_smooth[:,0]),len(oxy_m_I05_zave_smooth[0,:])+1)) 
oxy_m_I05_section[:,0] = oxy_m_I05_zave_smooth[:,0]
oxy_m_I05_section[:,2] = oxy_m_I05_zave_smooth[:,1]
for i,val in enumerate(oxy_m_I05_zave_smooth[:,0]):
    tmp = np.abs(val - lon_I05)
    oxy_m_I05_section[i,1] = lat_I05[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_I05_section = np.zeros((len(wma_m_I05_zave_smooth[:,0]),len(wma_m_I05_zave_smooth[0,:])+1)) 
wma_m_I05_section[:,0] = wma_m_I05_zave_smooth[:,0]
wma_m_I05_section[:,2] = wma_m_I05_zave_smooth[:,1]
for i,val in enumerate(wma_m_I05_zave_smooth[:,0]):
    tmp = np.abs(val - lon_I05)
    wma_m_I05_section[i,1] = lat_I05[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_I05_zave_nans = np.ma.getdata(oxy_o_I05_zave)
oxy_o_I05_zave_nans[oxy_o_I05_zave_nans<-10000] = np.nan
wma_o_I05_zave_nans = np.ma.getdata(wma_o_I05_zave)
wma_o_I05_zave_nans[wma_o_I05_zave_nans<-10000] = np.nan
oxy_o_I05_zave_smooth = lowess(oxy_o_I05_zave_nans, lon_I05, frac=0.2, delta=5)
wma_o_I05_zave_smooth = lowess(wma_o_I05_zave_nans, lon_I05, frac=0.2, delta=5)

oxy_o_I05_section = np.zeros((len(oxy_o_I05_zave_smooth[:,0]),len(oxy_o_I05_zave_smooth[0,:])+1)) 
oxy_o_I05_section[:,0] = oxy_o_I05_zave_smooth[:,0]
oxy_o_I05_section[:,2] = oxy_o_I05_zave_smooth[:,1]
for i,val in enumerate(oxy_o_I05_zave_smooth[:,0]):
    tmp = np.abs(val - lon_I05)
    oxy_o_I05_section[i,1] = lat_I05[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_I05_section = np.zeros((len(wma_o_I05_zave_smooth[:,0]),len(wma_o_I05_zave_smooth[0,:])+1)) 
wma_o_I05_section[:,0] = wma_o_I05_zave_smooth[:,0]
wma_o_I05_section[:,2] = wma_o_I05_zave_smooth[:,1]
for i,val in enumerate(wma_o_I05_zave_smooth[:,0]):
    tmp = np.abs(val - lon_I05)
    wma_o_I05_section[i,1] = lat_I05[int(np.where(tmp==np.min(tmp))[0][0])]


data.close()


#%% get data for A10

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_A10.nc', 'r')
oxy_m_A10 = np.squeeze(data.variables['AOU_M_DIF'][...])/22.5*10.
oxy_o_A10 = np.squeeze(data.variables['AOU_O_DIF'][...])/22.5*10.
wma_m_A10 = np.squeeze(data.variables['WMA_M_DIF'][...])/22.5*10.
wma_o_A10 = np.squeeze(data.variables['WMA_O_DIF'][...])/22.5*10.
oxy_m_A10_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/22.5*10.
oxy_o_A10_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/22.5*10.
wma_m_A10_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/22.5*10.
wma_o_A10_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/22.5*10.
rho_m_2011_A10 = np.squeeze(data.variables['RHO_M_2011'][...])
rho_o_2011_A10 = np.squeeze(data.variables['RHO_O_2011'][...])
rho_m_1988_1989_A10 = np.squeeze(data.variables['RHO_M_1988_1989'][...])
rho_o_1988_1989_A10 = np.squeeze(data.variables['RHO_O_1988_1989'][...])
oxy_m_rho_A10 = np.squeeze(data.variables['OXY_M_DIF_26P5_27P0'][...])/22.5*10.
oxy_o_rho_A10 = np.squeeze(data.variables['OXY_O_DIF_26P5_27P0'][...])/22.5*10.
sat_m_rho_A10 = np.squeeze(data.variables['SAT_M_DIF_26P5_27P0'][...])/22.5*10.
sat_o_rho_A10 = np.squeeze(data.variables['SAT_O_DIF_26P5_27P0'][...])/22.5*10.
aou_m_rho_A10 = np.squeeze(data.variables['AOU_M_DIF_26P5_27P0'][...])/22.5*10.*(-1)
aou_o_rho_A10 = np.squeeze(data.variables['AOU_O_DIF_26P5_27P0'][...])/22.5*10.*(-1)
age_m_rho_A10 = np.squeeze(data.variables['WMA_M_DIF_26P5_27P0'][...])/22.5*10.
age_o_rho_A10 = np.squeeze(data.variables['WMA_O_DIF_26P5_27P0'][...])/22.5*10.
dem_m_rho_A10 = np.squeeze(data.variables['DEM_M_DIF_26P5_27P0'][...])/22.5*10.
dem_o_rho_A10 = np.squeeze(data.variables['DEM_O_DIF_26P5_27P0'][...])/22.5*10.

sup_m_rho_A10 = aou_m_rho_A10 - dem_m_rho_A10
sup_o_rho_A10 = aou_o_rho_A10 - dem_o_rho_A10
fdem_m_rho_A10 = dem_m_rho_A10/aou_m_rho_A10
fsup_m_rho_A10 = sup_m_rho_A10/aou_m_rho_A10
fdem_o_rho_A10 = dem_o_rho_A10/aou_o_rho_A10
fsup_o_rho_A10 = sup_o_rho_A10/aou_o_rho_A10


dep_A10 = data.variables['DEPTHT11_23'][...]
lon_A10 = data.variables['ETOPO60X1N79_1'][...]
lat_A10 = np.zeros(len(lon_A10))
for i,val in enumerate(lon_A10):
    tmp = np.abs(val - a10[:,1])
    lat_A10[i] = a10[int(np.where(tmp==np.min(tmp))[0][0]),0]


# smoothing
oxy_m_A10_zave_nans = np.ma.getdata(oxy_m_A10_zave)
oxy_m_A10_zave_nans[oxy_m_A10_zave_nans<-10000] = np.nan
oxy_m_A10_zave_nans[oxy_m_A10_zave_nans<-50] = np.nan
wma_m_A10_zave_nans = np.ma.getdata(wma_m_A10_zave)
wma_m_A10_zave_nans[wma_m_A10_zave_nans<-10000] = np.nan
oxy_m_A10_zave_smooth = lowess(oxy_m_A10_zave_nans, lon_A10, frac=0.2, delta=5)
wma_m_A10_zave_smooth = lowess(wma_m_A10_zave_nans, lon_A10, frac=0.2, delta=5)

oxy_m_A10_section = np.zeros((len(oxy_m_A10_zave_smooth[:,0]),len(oxy_m_A10_zave_smooth[0,:])+1)) 
oxy_m_A10_section[:,0] = oxy_m_A10_zave_smooth[:,0]
oxy_m_A10_section[:,2] = oxy_m_A10_zave_smooth[:,1]
for i,val in enumerate(oxy_m_A10_zave_smooth[:,0]):
    tmp = np.abs(val - lon_A10)
    oxy_m_A10_section[i,1] = lat_A10[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_A10_section = np.zeros((len(wma_m_A10_zave_smooth[:,0]),len(wma_m_A10_zave_smooth[0,:])+1)) 
wma_m_A10_section[:,0] = wma_m_A10_zave_smooth[:,0]
wma_m_A10_section[:,2] = wma_m_A10_zave_smooth[:,1]
for i,val in enumerate(wma_m_A10_zave_smooth[:,0]):
    tmp = np.abs(val - lon_A10)
    wma_m_A10_section[i,1] = lat_A10[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_A10_zave_nans = np.ma.getdata(oxy_o_A10_zave)
oxy_o_A10_zave_nans[oxy_o_A10_zave_nans<-10000] = np.nan
wma_o_A10_zave_nans = np.ma.getdata(wma_o_A10_zave)
wma_o_A10_zave_nans[wma_o_A10_zave_nans<-10000] = np.nan
oxy_o_A10_zave_smooth = lowess(oxy_o_A10_zave_nans, lon_A10, frac=0.2, delta=5)
wma_o_A10_zave_smooth = lowess(wma_o_A10_zave_nans, lon_A10, frac=0.2, delta=5)

oxy_o_A10_section = np.zeros((len(oxy_o_A10_zave_smooth[:,0]),len(oxy_o_A10_zave_smooth[0,:])+1)) 
oxy_o_A10_section[:,0] = oxy_o_A10_zave_smooth[:,0]
oxy_o_A10_section[:,2] = oxy_o_A10_zave_smooth[:,1]
for i,val in enumerate(oxy_o_A10_zave_smooth[:,0]):
    tmp = np.abs(val - lon_A10)
    oxy_o_A10_section[i,1] = lat_A10[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_A10_section = np.zeros((len(wma_o_A10_zave_smooth[:,0]),len(wma_o_A10_zave_smooth[0,:])+1)) 
wma_o_A10_section[:,0] = wma_o_A10_zave_smooth[:,0]
wma_o_A10_section[:,2] = wma_o_A10_zave_smooth[:,1]
for i,val in enumerate(wma_o_A10_zave_smooth[:,0]):
    tmp = np.abs(val - lon_A10)
    wma_o_A10_section[i,1] = lat_A10[int(np.where(tmp==np.min(tmp))[0][0])]


data.close()



#%% get data for A13.5

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")
data = nc.Dataset('repeat_section_A13.5.nc', 'r')
oxy_m_A13 = np.squeeze(data.variables['AOU_M_DIF'][...])/26*10.
oxy_o_A13 = np.squeeze(data.variables['AOU_O_DIF'][...])/26*10.
wma_m_A13 = np.squeeze(data.variables['WMA_M_DIF'][...])/26*10.
wma_o_A13 = np.squeeze(data.variables['WMA_O_DIF'][...])/26*10.
oxy_m_A13_zave = np.squeeze(data.variables['AOU_M_DIF_ZAVE'][...])/26*10.
oxy_o_A13_zave = np.squeeze(data.variables['AOU_O_DIF_ZAVE'][...])/26*10.
wma_m_A13_zave = np.squeeze(data.variables['WMA_M_DIF_ZAVE'][...])/26*10.
wma_o_A13_zave = np.squeeze(data.variables['WMA_O_DIF_ZAVE'][...])/26*10.
rho_m_2010_A13 = np.squeeze(data.variables['RHO_M_2010'][...])
rho_o_2010_A13 = np.squeeze(data.variables['RHO_O_2010'][...])
rho_m_1983_1984_A13 = np.squeeze(data.variables['RHO_M_1983_1984'][...])
rho_o_1983_1984_A13 = np.squeeze(data.variables['RHO_O_1983_1984'][...])
oxy_m_rho_A13 = np.squeeze(data.variables['OXY_M_DIF_26P5_27P4'][...])/26*10.
oxy_o_rho_A13 = np.squeeze(data.variables['OXY_O_DIF_26P5_27P4'][...])/26*10.
sat_m_rho_A13 = np.squeeze(data.variables['SAT_M_DIF_26P5_27P4'][...])/26*10.
sat_o_rho_A13 = np.squeeze(data.variables['SAT_O_DIF_26P5_27P4'][...])/26*10.
aou_m_rho_A13 = np.squeeze(data.variables['AOU_M_DIF_26P5_27P4'][...])/26*10.*(-1)
aou_o_rho_A13 = np.squeeze(data.variables['AOU_O_DIF_26P5_27P4'][...])/26*10.*(-1)
age_m_rho_A13 = np.squeeze(data.variables['WMA_M_DIF_26P5_27P4'][...])/26*10.
age_o_rho_A13 = np.squeeze(data.variables['WMA_O_DIF_26P5_27P4'][...])/26*10.
dem_m_rho_A13 = np.squeeze(data.variables['DEM_M_DIF_26P5_27P4'][...])/26*10.
dem_o_rho_A13 = np.squeeze(data.variables['DEM_O_DIF_26P5_27P4'][...])/26*10.

sup_m_rho_A13 = aou_m_rho_A13 - dem_m_rho_A13
sup_o_rho_A13 = aou_o_rho_A13 - dem_o_rho_A13
fdem_m_rho_A13 = dem_m_rho_A13/aou_m_rho_A13
fsup_m_rho_A13 = sup_m_rho_A13/aou_m_rho_A13
fdem_o_rho_A13 = dem_o_rho_A13/aou_o_rho_A13
fsup_o_rho_A13 = sup_o_rho_A13/aou_o_rho_A13


dep_A13 = data.variables['DEPTHT11_23'][...]
lat_A13 = data.variables['ETOPO60Y121_100'][...]
lon_A13 = np.zeros(len(lat_A13))
for i,val in enumerate(lat_A13):
    tmp = np.abs(val - a135[:,0])
    lon_A13[i] = a135[int(np.where(tmp==np.min(tmp))[0][0]),1]



# smoothing
oxy_m_A13_zave_nans = np.ma.getdata(oxy_m_A13_zave)
oxy_m_A13_zave_nans[oxy_m_A13_zave_nans<-10000] = np.nan
wma_m_A13_zave_nans = np.ma.getdata(wma_m_A13_zave)
wma_m_A13_zave_nans[wma_m_A13_zave_nans<-10000] = np.nan
oxy_m_A13_zave_smooth = lowess(oxy_m_A13_zave_nans, lat_A13, frac=0.2, delta=5)
wma_m_A13_zave_smooth = lowess(wma_m_A13_zave_nans, lat_A13, frac=0.2, delta=5)

oxy_m_A13_section = np.zeros((len(oxy_m_A13_zave_smooth[:,0]),len(oxy_m_A13_zave_smooth[0,:])+1)) 
oxy_m_A13_section[:,1] = oxy_m_A13_zave_smooth[:,0]
oxy_m_A13_section[:,2] = oxy_m_A13_zave_smooth[:,1]
for i,val in enumerate(oxy_m_A13_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A13)
    oxy_m_A13_section[i,0] = lon_A13[int(np.where(tmp==np.min(tmp))[0][0])]

wma_m_A13_section = np.zeros((len(wma_m_A13_zave_smooth[:,0]),len(wma_m_A13_zave_smooth[0,:])+1)) 
wma_m_A13_section[:,1] = wma_m_A13_zave_smooth[:,0]
wma_m_A13_section[:,2] = wma_m_A13_zave_smooth[:,1]
for i,val in enumerate(wma_m_A13_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A13)
    wma_m_A13_section[i,0] = lon_A13[int(np.where(tmp==np.min(tmp))[0][0])]


oxy_o_A13_zave_nans = np.ma.getdata(oxy_o_A13_zave)
oxy_o_A13_zave_nans[oxy_o_A13_zave_nans<-10000] = np.nan
wma_o_A13_zave_nans = np.ma.getdata(wma_o_A13_zave)
wma_o_A13_zave_nans[wma_o_A13_zave_nans<-10000] = np.nan
oxy_o_A13_zave_smooth = lowess(oxy_o_A13_zave_nans, lat_A13, frac=0.2, delta=5)
wma_o_A13_zave_smooth = lowess(wma_o_A13_zave_nans, lat_A13, frac=0.2, delta=5)

oxy_o_A13_section = np.zeros((len(oxy_o_A13_zave_smooth[:,0]),len(oxy_o_A13_zave_smooth[0,:])+1)) 
oxy_o_A13_section[:,1] = oxy_o_A13_zave_smooth[:,0]
oxy_o_A13_section[:,2] = oxy_o_A13_zave_smooth[:,1]
for i,val in enumerate(oxy_o_A13_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A13)
    oxy_o_A13_section[i,0] = lon_A13[int(np.where(tmp==np.min(tmp))[0][0])]

wma_o_A13_section = np.zeros((len(wma_o_A13_zave_smooth[:,0]),len(wma_o_A13_zave_smooth[0,:])+1)) 
wma_o_A13_section[:,1] = wma_o_A13_zave_smooth[:,0]
wma_o_A13_section[:,2] = wma_o_A13_zave_smooth[:,1]
for i,val in enumerate(wma_o_A13_zave_smooth[:,0]):
    tmp = np.abs(val - lat_A13)
    wma_o_A13_section[i,0] = lon_A13[int(np.where(tmp==np.min(tmp))[0][0])]



data.close()


#%% figure specifics 

import matplotlib as mpl

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap = cmocean.tools.lighten(cmo.balance, 0.8)

levs1 = np.arange(-5,6,1)
levs2 = np.arange(-10,11,2)
norm1 = mpl.colors.BoundaryNorm(levs1, colmap.N)
norm2 = mpl.colors.BoundaryNorm(levs2, colmap.N)

contcol = 'black'
contwid = 0.75
alf = 0.8


lat_labs1 = ['80$^{\circ}$S', '60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq ']
lat_labs2 = [' ', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N']


#%%

def plotstations(p18, p16n, a135, a10, a16s, a16n, i05, cols,lws,alfs):
    #plt.plot(p18[:,1], p18[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)
    #plt.plot(p16n[:,1], p16n[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)
    plt.plot(a135[:,1], a135[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)
    plt.plot(a10[:,1], a10[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)
    plt.plot(a16s[:,1], a16s[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)
    plt.plot(a16n[:,1], a16n[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)
    plt.plot(i05[:,1], i05[:,0], transform=ccrs.PlateCarree(), color=cols, linewidth=lws, alpha=alfs)

cols = 'k'
lws = 1.5
alfs = 0.75



#%% make supp figure of IPSL-CM6A-LR trends

proj = ccrs.Robinson(central_longitude=30)
lons,lats = np.meshgrid(lon,lat)


fig = plt.figure(figsize=(12,4.5))
gs = GridSpec(2,2)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, aou_ipsl_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
plt.contour(lons, lats, aou_ipsl_tre, transform=ccrs.PlateCarree(), levels=np.array([0]), zorder=1, linewidths=lw, colors=contcol, alpha=alf)
plotstations(p18, p16n, a135, a10, a16s, a16n, i05, cols,lws,alfs)
ax1.set_extent([-90,120,-80,10])

ax2 = plt.subplot(gs[0,1], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, wma_ipsl_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='both')
plt.contour(lons, lats, wma_ipsl_tre, transform=ccrs.PlateCarree(), levels=np.array([0]), zorder=1, linewidths=lw, colors=contcol, alpha=alf)
plotstations(p18, p16n, a135, a10, a16s, a16n, i05, cols,lws,alfs)
ax2.set_extent([-90,120,-80,10])


ax3 = plt.subplot(gs[1,0], projection=proj)
ax3.tick_params(labelsize=fstic)
ax3.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax3.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax3.coastlines()
p3 = plt.scatter(oxy_o_I05_section[:,0], oxy_o_I05_section[:,1], c=oxy_o_I05_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, alpha=alf, norm=norm1)
plt.scatter(oxy_o_A13_section[:,0], oxy_o_A13_section[:,1], c=oxy_o_A13_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, alpha=alf, norm=norm1)
plt.scatter(oxy_o_A10_section[:,0], oxy_o_A10_section[:,1], c=oxy_o_A10_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, alpha=alf, norm=norm1)
plt.scatter(oxy_o_A16N_section[:,0], oxy_o_A16N_section[:,1], c=oxy_o_A16N_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, alpha=alf, norm=norm1)
plt.scatter(oxy_o_A16S_section[:,0], oxy_o_A16S_section[:,1], c=oxy_o_A16S_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, alpha=alf, norm=norm1)
ax3.set_extent([-90,120,-80,10])

ax4 = plt.subplot(gs[1,1], projection=proj)
ax4.tick_params(labelsize=fstic)
ax4.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax4.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax4.coastlines()
p4 = plt.scatter(wma_o_I05_section[:,0], wma_o_I05_section[:,1], c=wma_o_I05_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, alpha=alf, norm=norm2)
plt.scatter(wma_o_A13_section[:,0], wma_o_A13_section[:,1], c=wma_o_A13_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, alpha=alf, norm=norm2)
plt.scatter(wma_o_A10_section[:,0], wma_o_A10_section[:,1], c=wma_o_A10_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, alpha=alf, norm=norm2)
plt.scatter(wma_o_A16N_section[:,0], wma_o_A16N_section[:,1], c=wma_o_A16N_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, alpha=alf, norm=norm2)
plt.scatter(wma_o_A16S_section[:,0], wma_o_A16S_section[:,1], c=wma_o_A16S_section[:,2], transform=ccrs.PlateCarree(), cmap=colmap, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, alpha=alf, norm=norm2)
ax4.set_extent([-90,120,-80,10])


fig.subplots_adjust(top=0.95, bottom=0.1, left=0.125, right=0.875, wspace=0.1, hspace=0.1)


xx = 0.025; yy = 1.1
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'd', fontsize=fslab+2, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.09, 0.2, 0.025, 0.6])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=levs1[::2])
cbax1.set_ylabel('$\Delta$ AOU ($\mu$M decade$^{-1}$)', fontsize=fstic)
cbax1.tick_params(labelsize=fstic, right=False, labelright=False, left=True, labelleft=True)
cbax1.yaxis.set_label_position('left')

cbax2 = fig.add_axes([0.885, 0.2, 0.025, 0.6])
cbar2 = plt.colorbar(p2, cax=cbax2, orientation='vertical', ticks=levs2[::2])
cbax2.set_ylabel('$\Delta$ age (years decade$^{-1}$)', fontsize=fstic)
cbax2.tick_params(labelsize=fstic)


#%% save figure


os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig6.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig6_trans.png', dpi=300, bbox_inches='tight', transparent=True)

