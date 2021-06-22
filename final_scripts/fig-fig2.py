# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 10:39:44 2021

Purpose
-------
    Create figure of change in the b-value of the Martin curve, which has been 
    derived from the vertical profile of POC


@author: pearseb
"""

#%% imports

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sb
sb.set(style='ticks')

import cmocean
import cmocean.cm as cmo

from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy
import cartopy.crs as ccrs


#%% get data

os.chdir("C://Users/pearseb/Dropbox/PostDoc/my articles/historical model-data deoxygenation/data_for_figures")

data = nc.Dataset("figure_ageremin_vs_tou_z200-1000_ndep.nc",'r')
tou_tre = data.variables['TOU_TRE_ZAVE'][...]
wma_tre = data.variables['WMA_TRE_ZAVE'][...]
rem_tre = data.variables['RDEM_TRE_ZAVE'][...]
wma_r = data.variables['WMA_R'][...]
rem_r = data.variables['RDEM_R'][...]
wma_s = data.variables['WMA_SLOPE'][...]
rem_s = data.variables['RDEM_SLOPE'][...]
wma_std = data.variables['WMA_STD'][...]
rem_std = data.variables['RDEM_STD'][...]

deltaAOU_wma = wma_s * wma_tre
deltaAOU_rem = rem_s * rem_tre
deltaAOU_rem2wma = np.ma.abs(deltaAOU_rem) - np.ma.abs(deltaAOU_wma)

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()


#%% CMIP6 historical simulations

os.chdir("C://Users/pearseb/Dropbox/PostDoc/my articles/historical model-data deoxygenation/data_for_figures")

data = nc.Dataset("ACCESS-ESM1-5_historical_trends_and_correlations.nc",'r')
wma_r_access = data.variables['WMA_R'][...]
npp_r_access = data.variables['NPP_R'][...]
wma_s_access = data.variables['WMA_SLOPE'][...]
npp_s_access = data.variables['NPP_SLOPE'][...]
wma_tre_access = data.variables['WMA_TRE'][...]
npp_tre_access = data.variables['NPP_TRE'][...]

data = nc.Dataset("CanESM5_historical_trends_and_correlations.nc",'r')
wma_r_canesm = data.variables['WMA_R'][...]
npp_r_canesm = data.variables['NPP_R'][...]
wma_s_canesm = data.variables['WMA_SLOPE'][...]
npp_s_canesm = data.variables['NPP_SLOPE'][...]
wma_tre_canesm = data.variables['WMA_TRE'][...]
npp_tre_canesm = data.variables['NPP_TRE'][...]

data = nc.Dataset("IPSL-CM6A-LR_historical_trends_and_correlations.nc",'r')
wma_r_ipsl = data.variables['WMA_R'][...]
npp_r_ipsl = data.variables['NPP_R'][...]
wma_s_ipsl = data.variables['WMA_SLOPE'][...]
npp_s_ipsl = data.variables['NPP_SLOPE'][...]
wma_tre_ipsl = data.variables['WMA_TRE'][...]
npp_tre_ipsl = data.variables['NPP_TRE'][...]

data = nc.Dataset("MIROC-ES2L_historical_trends_and_correlations.nc",'r')
wma_r_miroc = data.variables['WMA_R'][...]
npp_r_miroc = data.variables['NPP_R'][...]
wma_s_miroc = data.variables['WMA_SLOPE'][...]
npp_s_miroc = data.variables['NPP_SLOPE'][...]
wma_tre_miroc = data.variables['WMA_TRE'][...]
npp_tre_miroc = data.variables['NPP_TRE'][...]

data = nc.Dataset("MPI-ESM1-2-LR_historical_trends_and_correlations.nc",'r')
wma_r_mpi = data.variables['WMA_R'][...]
npp_r_mpi = data.variables['NPP_R'][...]
wma_s_mpi = data.variables['WMA_SLOPE'][...]
npp_s_mpi = data.variables['NPP_SLOPE'][...]
wma_tre_mpi = data.variables['WMA_TRE'][...]
npp_tre_mpi = data.variables['NPP_TRE'][...]

data = nc.Dataset("MRI-ESM2-0_historical_trends_and_correlations.nc",'r')
wma_r_mri = data.variables['WMA_R'][...]
npp_r_mri = data.variables['NPP_R'][...]
wma_s_mri = data.variables['WMA_SLOPE'][...]
npp_s_mri = data.variables['NPP_SLOPE'][...]
wma_tre_mri = data.variables['WMA_TRE'][...]
npp_tre_mri = data.variables['NPP_TRE'][...]

data = nc.Dataset("UKESM1-0-LL_historical_trends_and_correlations.nc",'r')
wma_r_ukesm = data.variables['WMA_R'][...]
npp_r_ukesm = data.variables['NPP_R'][...]
wma_s_ukesm = data.variables['WMA_SLOPE'][...]
npp_s_ukesm = data.variables['NPP_SLOPE'][...]
wma_tre_ukesm = data.variables['WMA_TRE'][...]
npp_tre_ukesm = data.variables['NPP_TRE'][...]


#%% calculate multi-model mean correlations and agreement

wma_r_cmip6 = np.ma.zeros((7,180,360))
wma_r_cmip6[0,:,:] = wma_r_access
wma_r_cmip6[1,:,:] = wma_r_canesm
wma_r_cmip6[2,:,:] = wma_r_ipsl
wma_r_cmip6[3,:,:] = wma_r_miroc
wma_r_cmip6[4,:,:] = wma_r_mpi
wma_r_cmip6[5,:,:] = wma_r_mri
wma_r_cmip6[6,:,:] = wma_r_ukesm

npp_r_cmip6 = np.ma.zeros((7,180,360))
npp_r_cmip6[0,:,:] = npp_r_access
npp_r_cmip6[1,:,:] = npp_r_canesm
npp_r_cmip6[2,:,:] = npp_r_ipsl
npp_r_cmip6[3,:,:] = npp_r_miroc
npp_r_cmip6[4,:,:] = npp_r_mpi
npp_r_cmip6[5,:,:] = npp_r_mri
npp_r_cmip6[6,:,:] = npp_r_ukesm

wma_r_cmip6_mean = np.ma.average(wma_r_cmip6, axis=0)
npp_r_cmip6_mean = np.ma.average(npp_r_cmip6, axis=0)

# calculate if there is at least 80% agreement between models on sign | at least 6 out of the 7 models must agree
wma_r_cmip6_sign = np.sign(wma_r_cmip6)
wma_r_cmip6_coun = np.ma.count(wma_r_cmip6, axis=0)
wma_r_cmip6_pos = np.ma.masked_where(wma_r_cmip6_sign < 0, wma_r_cmip6_sign)
wma_r_cmip6_neg = np.ma.masked_where(wma_r_cmip6_sign > 0, wma_r_cmip6_sign)
wma_r_cmip6_pfra = np.ma.sum(wma_r_cmip6_pos, axis=0) / wma_r_cmip6_coun
wma_r_cmip6_nfra = np.ma.sum(wma_r_cmip6_neg, axis=0) / wma_r_cmip6_coun 
wma_r_cmip6_agre = np.ma.zeros((180,360))
wma_r_cmip6_agre[wma_r_cmip6_pfra > 0.8] = 1
wma_r_cmip6_agre[wma_r_cmip6_nfra < -0.8] = 1

# calculate if there is at least 80% agreement between models on sign | at least 6 out of the 7 models must agree
npp_r_cmip6_sign = np.sign(npp_r_cmip6)
npp_r_cmip6_coun = np.ma.count(npp_r_cmip6, axis=0)
npp_r_cmip6_pos = np.ma.masked_where(npp_r_cmip6_sign < 0, npp_r_cmip6_sign)
npp_r_cmip6_neg = np.ma.masked_where(npp_r_cmip6_sign > 0, npp_r_cmip6_sign)
npp_r_cmip6_pfra = np.ma.sum(npp_r_cmip6_pos, axis=0) / npp_r_cmip6_coun
npp_r_cmip6_nfra = np.ma.sum(npp_r_cmip6_neg, axis=0) / npp_r_cmip6_coun 
npp_r_cmip6_agre = np.ma.zeros((180,360))
npp_r_cmip6_agre[npp_r_cmip6_pfra > 0.8] = 1
npp_r_cmip6_agre[npp_r_cmip6_nfra < -0.8] = 1


#%% calculate the difference between Age- and NPP-driven AOU changes in the mesopelagic

deltaAOU_wma_access = wma_s_access * wma_tre_access
deltaAOU_wma_canesm = wma_s_canesm * wma_tre_canesm
deltaAOU_wma_ipsl = wma_s_ipsl * wma_tre_ipsl
deltaAOU_wma_miroc = wma_s_miroc * wma_tre_miroc
deltaAOU_wma_mpi = wma_s_mpi * wma_tre_mpi
deltaAOU_wma_mri = wma_s_mri * wma_tre_mri
deltaAOU_wma_ukesm = wma_s_ukesm * wma_tre_ukesm

deltaAOU_npp_access = npp_s_access * npp_tre_access
deltaAOU_npp_canesm = npp_s_canesm * npp_tre_canesm
deltaAOU_npp_ipsl = npp_s_ipsl * npp_tre_ipsl
deltaAOU_npp_miroc = npp_s_miroc * npp_tre_miroc
deltaAOU_npp_mpi = npp_s_mpi * npp_tre_mpi
deltaAOU_npp_mri = npp_s_mri * npp_tre_mri
deltaAOU_npp_ukesm = npp_s_ukesm * npp_tre_ukesm


deltaAOU_wma_cmip6 = np.ma.zeros((7,180,360))
deltaAOU_wma_cmip6[0,:,:] = deltaAOU_wma_access
deltaAOU_wma_cmip6[1,:,:] = deltaAOU_wma_canesm
deltaAOU_wma_cmip6[2,:,:] = deltaAOU_wma_ipsl
deltaAOU_wma_cmip6[3,:,:] = deltaAOU_wma_miroc
deltaAOU_wma_cmip6[4,:,:] = deltaAOU_wma_mpi
deltaAOU_wma_cmip6[5,:,:] = deltaAOU_wma_mri
deltaAOU_wma_cmip6[6,:,:] = deltaAOU_wma_ukesm

deltaAOU_npp_cmip6 = np.ma.zeros((7,180,360))
deltaAOU_npp_cmip6[0,:,:] = deltaAOU_npp_access
deltaAOU_npp_cmip6[1,:,:] = deltaAOU_npp_canesm
deltaAOU_npp_cmip6[2,:,:] = deltaAOU_npp_ipsl
deltaAOU_npp_cmip6[3,:,:] = deltaAOU_npp_miroc
deltaAOU_npp_cmip6[4,:,:] = deltaAOU_npp_mpi
deltaAOU_npp_cmip6[5,:,:] = deltaAOU_npp_mri
deltaAOU_npp_cmip6[6,:,:] = deltaAOU_npp_ukesm


deltaAOU_npp2wma_cmip6 = np.ma.abs(deltaAOU_npp_cmip6) - np.ma.abs(deltaAOU_wma_cmip6)
deltaAOU_npp2wma_cmip6_mean = np.ma.average(deltaAOU_npp2wma_cmip6, axis=0)


# calculate the standard deviation, and see if it crosses zero
deltaAOU_npp2wma_cmip6_1std = np.ma.std(deltaAOU_npp2wma_cmip6, axis=0)
deltaAOU_npp2wma_cmip6_agre = np.ma.zeros((180,360))
deltaAOU_npp2wma_cmip6_agre[np.ma.abs(deltaAOU_npp2wma_cmip6_mean) > deltaAOU_npp2wma_cmip6_1std] = 1


#%% figure specifics

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.balance, 0.8)
levs1 = np.array([-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5])
colmap2 = cmocean.tools.lighten(cmo.balance, 0.8)
levs2 = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

cont1 = np.array([-2, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
cont1 = np.array([-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5])


contcol = 'black'
contwid=0.75

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)



#%% figure

fig = plt.figure(figsize=(14,6))
gs = GridSpec(2,2)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, wma_r, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='neither')
c1 = plt.contour(lons, lats, wma_r, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax2 = plt.subplot(gs[0,1], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, rem_r, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='neither')
c2 = plt.contour(lons, lats, rem_r, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

'''
ax3 = plt.subplot(gs[0,2], projection=proj)
ax3.tick_params(labelsize=fstic)
ax3.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax3.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax3.coastlines()
p3 = plt.contourf(lons, lats, deltaAOU_rem2wma, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='neither')
#c3 = plt.contour(lons, lats, deltaAOU_rem2wma, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)
#plt.contourf(lons, lats, rem_r, transform=ccrs.PlateCarree(), colors='none', levels=[-1,0.2,1], hatches=['............',' '], zorder=2)
'''

ax4 = plt.subplot(gs[1,0], projection=proj)
ax4.tick_params(labelsize=fstic)
ax4.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax4.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax4.coastlines()
p4 = plt.contourf(lons, lats, wma_r_cmip6_mean, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='neither')
c4 = plt.contour(lons, lats, wma_r_cmip6_mean, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)
plt.contourf(lons, lats, wma_r_cmip6_agre, transform=ccrs.PlateCarree(), colors='none', levels=[0.5,1.5], hatches=['.....'], zorder=3)


ax5 = plt.subplot(gs[1,1], projection=proj)
ax5.tick_params(labelsize=fstic)
ax5.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax5.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax5.coastlines()
p5 = plt.contourf(lons, lats, npp_r_cmip6_mean, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='both')
c5 = plt.contour(lons, lats, npp_r_cmip6_mean, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)
plt.contourf(lons, lats, npp_r_cmip6_agre, transform=ccrs.PlateCarree(), colors='none', levels=[0.5,1.5], hatches=['.....'], zorder=3)

'''
ax6 = plt.subplot(gs[1,2], projection=proj)
ax6.tick_params(labelsize=fstic)
ax6.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax6.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax6.coastlines()
p6 = plt.contourf(lons, lats, deltaAOU_npp2wma_cmip6_mean, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='both')
#c6 = plt.contour(lons, lats, deltaAOU_npp2wma_cmip6_mean, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)
#plt.contourf(lons, lats, npp_r_cmip6_mean, transform=ccrs.PlateCarree(), colors='none', levels=[-1,0.2,1], hatches=['............',' '], zorder=2)
'''

fig.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9, wspace=-0.1, hspace=0.15)



xx = -0.1; yy = 0.5
plt.text(xx, yy, 'Hindcast simulation', fontsize=fslab, transform=ax1.transAxes, va='center', ha='center', fontweight='bold', rotation=90)
plt.text(xx, yy, 'CMIP6 Historical', fontsize=fslab, transform=ax4.transAxes, va='center', ha='center', fontweight='bold', rotation=90)

xx = 0.5; yy = 1.05
plt.text(xx, yy, 'Ideal age', fontsize=fstic, transform=ax1.transAxes, va='center', ha='center')
plt.text(xx, yy, 'Ideal age', fontsize=fstic, transform=ax4.transAxes, va='center', ha='center')
plt.text(xx, yy, 'Biological demand', fontsize=fstic, transform=ax2.transAxes, va='center', ha='center')
plt.text(xx, yy, 'Net Primary Production', fontsize=fstic, transform=ax5.transAxes, va='center', ha='center')
#plt.text(xx, yy, 'effect of demand relative to age', fontsize=fstic, transform=ax3.transAxes, va='center', ha='center')
#plt.text(xx, yy, 'effect of NPP relative to age', fontsize=fstic, transform=ax6.transAxes, va='center', ha='center')

xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
#plt.text(xx, yy, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'c', fontsize=fslab+2, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'd', fontsize=fslab+2, transform=ax5.transAxes, va='center', ha='center', fontweight='bold')
#plt.text(xx, yy, 'f', fontsize=fslab+2, transform=ax6.transAxes, va='center', ha='center', fontweight='bold')



cbax1 = fig.add_axes([0.9, 0.2, 0.03, 0.6])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=levs2[::2])
cbax1.tick_params(labelsize=fstic, bottom=False, labelbottom=False, top=True, labeltop=True)
cbax1.set_ylabel("Pearson's correlation with AOU", fontsize=fslab)
cbax1.yaxis.set_label_position('right')

'''
cbax2 = fig.add_axes([0.435, 0.85, 0.18, 0.035])
cbar2 = plt.colorbar(p3, cax=cbax2, orientation='horizontal', ticks=levs2[::2])
cbax2.tick_params(labelsize=fstic, bottom=False, labelbottom=False, top=True, labeltop=True)
cbax2.set_xlabel("correlation with AOU", fontsize=fslab)
cbax2.xaxis.set_label_position('top')

cbax3 = fig.add_axes([0.725, 0.85, 0.18, 0.035])
cbar3 = plt.colorbar(p5, cax=cbax3, orientation='horizontal', ticks=levs2[::2])
cbax3.tick_params(labelsize=fstic, bottom=False, labelbottom=False, top=True, labeltop=True)
cbax3.set_xlabel("$\mu$M AOU decade$^{-1}$", fontsize=fslab)
cbax3.xaxis.set_label_position('top')
'''



#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-fig2.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-fig2_trans.png', dpi=300, bbox_inches='tight', transparent=True)


