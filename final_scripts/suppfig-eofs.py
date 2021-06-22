# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 10:08:25 2021

Purpose
-------
    EOFs of oxygen and preformed oxygen output from hindcast simulation


@author: pearseb
"""

#%% imports

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
import cmocean
import cmocean.cm as cmo
import mpl_toolkits.basemap as bm
import datetime

import seaborn as sb
sb.set(style='ticks')


#%% get hindcast simulation

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")

data = nc.Dataset('ETOPO_JRA55_ndep_1m_oxygen4EOFanalysis_T_1972-2019.nc', 'r')
print(data)
o2_meso = np.ma.squeeze(data.variables['O2_MESO'][...])
year = np.arange(1972.5,2020,1)
mnth = np.arange(1,13,1)

lat = data.variables['ETOPO60Y'][...]
lon = data.variables['ETOPO60X'][...]
lons_np, lats_np = np.meshgrid(lon,lat)


#%% remove the seasonal cycle from the data

o2_meso = o2_meso - np.ma.average(o2_meso, axis=0)


#%% put all variables onto monthly timeseries and mask everywhere not in North Pacific

### put all data on on only one time axis
o2_meso_mnth = np.ma.zeros((48*12,180,360))
for ii,yr in enumerate(year):
    for jj in np.arange(12):
        idx = ii*12+jj
        o2_meso_mnth[idx] = o2_meso[ii,jj,:,:]


lons3D = np.repeat(lons_np[np.newaxis, :, :], 12*48, axis=0)
lats3D = np.repeat(lats_np[np.newaxis, :, :], 12*48, axis=0)

domain_np = [20,130,65,260]                

# isolate the North Pacific region
o2_meso_np = np.ma.masked_where(lons3D<domain_np[1],o2_meso_mnth)
o2_meso_np = np.ma.masked_where(lons3D>domain_np[3],o2_meso_np)
o2_meso_np = np.ma.masked_where(lats3D<domain_np[0],o2_meso_np)
o2_meso_np = np.ma.masked_where(lats3D>domain_np[2],o2_meso_np)


#%% North Atlantic and Southern Ocean

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")

data = nc.Dataset('ETOPO_JRA55_ndep_1m_oxygen4EOFanalysis_T_1972-2019.nc', 'r')
print(data)
o2_meso = np.ma.squeeze(data.variables['O2_MESO'][...])
year = np.arange(1972.5,2020,1)
mnth = np.arange(1,13,1)
lat = data.variables['ETOPO60Y'][...]
lon = data.variables['ETOPO60X'][...]

lon[lon>360] -= 360.0
lon = np.ma.concatenate((lon[340:360], lon[0:340]), axis=0)
o2_meso = np.ma.concatenate((o2_meso[:,:,:,340:360], o2_meso[:,:,:,0:340]), axis=3)

lon[lon>180.0] -= 360.0
lon = np.ma.concatenate((lon[180:360], lon[0:180]), axis=0)
o2_meso = np.ma.concatenate((o2_meso[:,:,:,180:360], o2_meso[:,:,:,0:180]), axis=3)

lons,lats = np.meshgrid(lon,lat)
lons = np.ma.concatenate((np.reshape(lons[:,359],(180,1)), lons[:,:]), axis=1)
lons[:,0] -= 360.0
lats = np.ma.concatenate((np.reshape(lats[:,-1], (180,1)), lats[:,:]), axis=1)
o2_meso = np.ma.concatenate((np.reshape(o2_meso[:,:,:,-1], (48,12,180,1)), o2_meso[:,:,:,:]), axis=3)

data.close()


#%% remove the seasonal cycle from the data

o2_meso = o2_meso - np.ma.average(o2_meso, axis=0)


#%% put all variables onto monthly timeseries and mask everywhere not in North Pacific

### put all data on on only one time axis
o2_meso_mnth = np.ma.zeros((48*12,180,361))
for ii,yr in enumerate(year):
    for jj in np.arange(12):
        idx = ii*12+jj
        o2_meso_mnth[idx] = o2_meso[ii,jj,:,:]


lons3D = np.repeat(lons[np.newaxis, :, :], 12*48, axis=0)
lats3D = np.repeat(lats[np.newaxis, :, :], 12*48, axis=0)

domain_so = [-80,-180,0,180]                

# isolate the Southern Ocean region
o2_meso_so = np.ma.masked_where(lons3D<domain_so[1],o2_meso_mnth)
o2_meso_so = np.ma.masked_where(lons3D>domain_so[3],o2_meso_so)
o2_meso_so = np.ma.masked_where(lats3D<domain_so[0],o2_meso_so)
o2_meso_so = np.ma.masked_where(lats3D>domain_so[2],o2_meso_so)

domain_na = [20,-90,80,20]                

# isolate the North Atlantic region
o2_meso_na = np.ma.masked_where(lons3D<domain_na[1],o2_meso_mnth)
o2_meso_na = np.ma.masked_where(lons3D>domain_na[3],o2_meso_na)
o2_meso_na = np.ma.masked_where(lats3D<domain_na[0],o2_meso_na)
o2_meso_na = np.ma.masked_where(lats3D>domain_na[2],o2_meso_na)



#%% EOF analysis

from eofs.standard import Eof

# create weights
coslat = np.cos(np.deg2rad(lats3D[:,:,0:360])).clip(0,1.0)
wgts = np.sqrt(coslat)

# solve EOFs
solver_o2_np = Eof(o2_meso_np[:,:,:], weights=wgts)
eofs_o2_np = solver_o2_np.eofsAsCovariance(neofs=3)
pcs_o2_np = solver_o2_np.pcs(npcs=3, pcscaling=1)
eigens_o2_np = solver_o2_np.eigenvalues(neigs=3)
varfrac_o2_np = solver_o2_np.varianceFraction(neigs=3)

coslat = np.cos(np.deg2rad(lats3D[:,:,:])).clip(0,1.0)
wgts = np.sqrt(coslat)

solver_o2_na = Eof(o2_meso_na[:,:,:], weights=wgts)
eofs_o2_na = solver_o2_na.eofsAsCovariance(neofs=3)
pcs_o2_na = solver_o2_na.pcs(npcs=3, pcscaling=1)
eigens_o2_na = solver_o2_na.eigenvalues(neigs=3)
varfrac_o2_na = solver_o2_na.varianceFraction(neigs=3)

solver_o2_so = Eof(o2_meso_so[:,:,:], weights=wgts)
eofs_o2_so = solver_o2_so.eofsAsCovariance(neofs=3)
pcs_o2_so = solver_o2_so.pcs(npcs=3, pcscaling=1)
eigens_o2_so = solver_o2_so.eigenvalues(neigs=3)
varfrac_o2_so = solver_o2_so.varianceFraction(neigs=3)



#%% get the PDO, AMO, SAM index from the hindcast simulation

os.chdir("C://Users//pearseb//Dropbox//PostDoc//Data//climate_modes")
pdo = pd.read_csv('pdo_noaa.csv', header=1)
pdo[pdo == -9999] = np.nan
dates_pdo = pd.date_range(start="1854-01-15", end="2021-08-15", periods=2012)

os.chdir("C://Users//pearseb//Dropbox//PostDoc//Data//climate_modes")
amo = pd.read_csv('amo.csv', header=0)
amo[amo == -9999] = np.nan
dates_amo = pd.date_range(start="1958-01-15", end="2019-12-15", periods=744)

os.chdir("C://Users//pearseb//Dropbox//PostDoc//Data//climate_modes")
sam = pd.read_csv('sam.csv', header=0)
dates_sam = pd.date_range(start="1979-06-14", end="2019-06-14", periods=41)



#%% view EOFs of mesopelagic O2

levs1 = np.arange(-40,41,4)*0.1
levs2 = np.arange(-40,41,4)*0.1
colmap1 = cmocean.tools.lighten(cmo.balance, 0.75)
colmap2 = cmocean.tools.lighten(cmo.balance, 0.75)
fstic = 14
fslab = 16

fig = plt.figure(figsize=(12,9.5), facecolor='w')
gs = GridSpec(22,2)

domain = [15,120,70,270]                
domain_draw = [20,130,65,260]
dlat=10
dlon=30

proj = bm.Basemap(projection='mill', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
lonproj, latproj = proj(lons_np, lats_np)


ax1 = plt.subplot(gs[0:8,0])
ax1.tick_params(labelsize=fstic)
proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
p1 = plt.contourf(lonproj,latproj, eofs_o2_np[0,:,:], cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
proj.drawparallels(range(domain_draw[0], domain_draw[2]+1, dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
proj.drawmeridians(range(domain_draw[1], domain_draw[3]+1, dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
plt.colorbar(p1, orientation='horizontal')

ax2 = plt.subplot(gs[0:8,1])
ax2.tick_params(labelsize=fstic)
proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
p2 = plt.contourf(lonproj,latproj, eofs_o2_np[1,:,:], cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
proj.drawparallels(range(domain_draw[0], domain_draw[2]+1, dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
proj.drawmeridians(range(domain_draw[1], domain_draw[3]+1, dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
plt.colorbar(p2, orientation='horizontal')


dates = pd.date_range(start="1972-01-15", end="2019-12-15", periods=576)

ax3 = plt.subplot(gs[8:10,0])
ax3.tick_params(labelsize=fstic)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.plot(dates, pcs_o2_np[:,0])
pcs_pos = np.ma.masked_where(pcs_o2_np[:,0]<0, pcs_o2_np[:,0])
pcs_neg = np.ma.masked_where(pcs_o2_np[:,0]>0, pcs_o2_np[:,0])
plt.fill_between(dates, pcs_pos, alpha=0.5, facecolor='firebrick')
plt.fill_between(dates, pcs_neg, alpha=0.5, facecolor='royalblue')
plt.plot(dates_pdo[60::], pdo['JRA55'][0:-60]*(-1), color='k', linewidth=1.5, alpha=0.5, linestyle='-')
plt.fill_betweenx((-2.5,2.5), datetime.date(1975,1,1), datetime.date(1984,12,31), alpha=0.2, color='k', edgecolor='none')
plt.fill_betweenx((-2.5,2.5), datetime.date(2005,1,1), datetime.date(2014,12,31), alpha=0.2, color='k', edgecolor='none')
plt.xlim(dates[0],dates[-1])

ax4 = plt.subplot(gs[8:10,1])
ax4.tick_params(labelsize=fstic)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
plt.plot(dates, pcs_o2_np[:,1])
pcs_pos = np.ma.masked_where(pcs_o2_np[:,1]<0, pcs_o2_np[:,1])
pcs_neg = np.ma.masked_where(pcs_o2_np[:,1]>0, pcs_o2_np[:,1])
plt.fill_between(dates, pcs_pos, alpha=0.5, facecolor='firebrick')
plt.fill_between(dates, pcs_neg, alpha=0.5, facecolor='royalblue')
plt.plot(dates_pdo[60::], pdo['JRA55'][0:-60]*(-1), color='k', linewidth=1.5, alpha=0.5, linestyle='-')
plt.fill_betweenx((-2.5,2.5), datetime.date(1975,1,1), datetime.date(1984,12,31), alpha=0.2, color='k', edgecolor='none')
plt.fill_betweenx((-2.5,2.5), datetime.date(2005,1,1), datetime.date(2014,12,31), alpha=0.2, color='k', edgecolor='none')
plt.xlim(dates[0],dates[-1])




domain = [20,-90,80,20]                
domain_draw = [20,-90,80,20]
dlat=10
dlon=30

proj = bm.Basemap(projection='mill', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
lonproj, latproj = proj(lons, lats)

ax5 = plt.subplot(gs[12:20,0])
ax5.tick_params(labelsize=fstic)
proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
p5 = plt.contourf(lonproj,latproj, eofs_o2_na[0,:,:], cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
proj.drawparallels(range(domain_draw[0], domain_draw[2]+1, dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
proj.drawmeridians(range(domain_draw[1], domain_draw[3]+1, dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
plt.colorbar(p5, orientation='horizontal')

domain = [-80,-180,0,180]                
domain_draw = [-80,-180,0,180]
dlat=20
dlon=60

proj = bm.Basemap(projection='mill', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
lonproj, latproj = proj(lons, lats)

ax6 = plt.subplot(gs[12:20,1])
ax6.tick_params(labelsize=fstic)
proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
p6 = plt.contourf(lonproj,latproj, eofs_o2_so[0,:,:], cmap=colmap2, levels=levs2, vmin=np.min(levs2), vmax=np.max(levs2), extend='both')
proj.drawparallels(range(domain_draw[0], domain_draw[2]+1, dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
proj.drawmeridians(range(domain_draw[1], domain_draw[3]+1, dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=fstic)
plt.colorbar(p6, orientation='horizontal')

ax7 = plt.subplot(gs[20:22,0])
ax7.tick_params(labelsize=fstic)
ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
plt.plot(dates, pcs_o2_na[:,0])
pcs_pos = np.ma.masked_where(pcs_o2_na[:,0]<0, pcs_o2_na[:,0])
pcs_neg = np.ma.masked_where(pcs_o2_na[:,0]>0, pcs_o2_na[:,0])
plt.fill_between(dates, pcs_pos, alpha=0.5, facecolor='firebrick')
plt.fill_between(dates, pcs_neg, alpha=0.5, facecolor='royalblue')
plt.plot(dates_amo, amo['JRA55']*10, color='k', linewidth=1.5, alpha=0.5, linestyle='-')
plt.fill_betweenx((-4,4), datetime.date(1975,1,1), datetime.date(1984,12,31), alpha=0.2, color='k', edgecolor='none')
plt.fill_betweenx((-4,4), datetime.date(2005,1,1), datetime.date(2014,12,31), alpha=0.2, color='k', edgecolor='none')
plt.xlim(dates[0],dates[-1])

ax8 = plt.subplot(gs[20:22,1])
ax8.tick_params(labelsize=fstic)
ax8.spines['top'].set_visible(False)
ax8.spines['right'].set_visible(False)
plt.plot(dates, pcs_o2_so[:,0])
pcs_pos = np.ma.masked_where(pcs_o2_so[:,0]<0, pcs_o2_so[:,0])
pcs_neg = np.ma.masked_where(pcs_o2_so[:,0]>0, pcs_o2_so[:,0])
plt.fill_between(dates, pcs_pos, alpha=0.5, facecolor='firebrick')
plt.fill_between(dates, pcs_neg, alpha=0.5, facecolor='royalblue')
plt.plot(dates_sam, sam['JRA55'], color='k', linewidth=1.5, alpha=0.5, linestyle='-')
plt.fill_betweenx((-6,6), datetime.date(1975,1,1), datetime.date(1984,12,31), alpha=0.2, color='k', edgecolor='none')
plt.fill_betweenx((-6,6), datetime.date(2005,1,1), datetime.date(2014,12,31), alpha=0.2, color='k', edgecolor='none')
plt.xlim(dates[0],dates[-1])


plt.text(0,1.05,'EOF1 - %.2f%% of variability'%(varfrac_o2_np[0]*100), transform=ax1.transAxes)
plt.text(0,1.05,'EOF2 - %.2f%% of variability'%(varfrac_o2_np[1]*100), transform=ax2.transAxes)
plt.text(0,1.05,'EOF1 - %.2f%% of variability'%(varfrac_o2_na[0]*100), transform=ax5.transAxes)
plt.text(0,1.05,'EOF1 - %.2f%% of variability'%(varfrac_o2_so[0]*100), transform=ax6.transAxes)


ax3.set_ylim(-2.5,2.5)
ax4.set_ylim(-2.5,2.5)
ax7.set_ylim(-3.5,3.5)
ax8.set_ylim(-5,5)

xx = -0.075; yy = 1.0
plt.text(xx,yy, 'a', transform=ax1.transAxes, va='center', ha='center', fontsize=fslab+2, fontweight='bold')
plt.text(xx,yy, 'b', transform=ax2.transAxes, va='center', ha='center', fontsize=fslab+2, fontweight='bold')
plt.text(xx,yy, 'c', transform=ax5.transAxes, va='center', ha='center', fontsize=fslab+2, fontweight='bold')
plt.text(xx,yy, 'd', transform=ax6.transAxes, va='center', ha='center', fontsize=fslab+2, fontweight='bold')


#%% save figure

os.chdir("C://Users//pearseb//Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('suppfig-eofs.png', dpi=300, bbox_inches='tight')
fig.savefig('suppfig-eofs_trans.png', dpi=300, bbox_inches='tight', transparent=True)
