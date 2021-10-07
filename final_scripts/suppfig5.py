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

data = nc.Dataset("ETOPO_budget_supplyVdemand.nc",'r')
sup = data.variables['DELTA_SUP'][...]*365
dem = data.variables['DELTA_DEM'][...]*365
suppro = data.variables['DELTA_SUP_PRO'][...]

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()



#%% figure specifics

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.balance, 0.8)
levs1 = np.arange(-1, 1.1, 0.1)*1e-1

colmap2 = cmocean.tools.lighten(cmo.deep, 0.8)
levs2 = np.arange(0.5, 1.01, 0.05)

contcol = 'black'
contwid=0.75

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)



#%% figure

fig = plt.figure(figsize=(7,9.5))
gs = GridSpec(3,1)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, sup, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
c1 = plt.contour(lons, lats, sup, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax2 = plt.subplot(gs[1,0], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, dem, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
c2 = plt.contour(lons, lats, dem, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax3 = plt.subplot(gs[2,0], projection=proj)
ax3.tick_params(labelsize=fstic)
ax3.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax3.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax3.coastlines()
p3 = plt.contourf(lons, lats, suppro, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=2, extend='min')
c3 = plt.contour(lons, lats, suppro, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

fig.subplots_adjust(top=0.95, bottom=0.05, left=0.05, right=0.82, wspace=0.05, hspace=0.2)


xx = 0.5; yy = 1.05
plt.text(xx, yy, 'Supply (advection + mixing)', fontsize=fslab, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'Demand (aerobic processes)', fontsize=fslab, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'O$_2$ flux due to supply', fontsize=fslab, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')

xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.82, 0.6925, 0.03, 0.25])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=levs1[::4])
cbax1.tick_params(labelsize=fstic)
cbax1.set_ylabel("$\Delta$ O$_2$ ($\mu$M year$^{-1}$)", fontsize=fslab)

cbax2 = fig.add_axes([0.82, 0.375, 0.03, 0.25])
cbar2 = plt.colorbar(p2, cax=cbax2, orientation='vertical', ticks=levs1[::4])
cbax2.tick_params(labelsize=fstic)
cbax2.set_ylabel("$\Delta$ O$_2$ ($\mu$M year$^{-1}$)", fontsize=fslab)

cbax3 = fig.add_axes([0.82, 0.0575, 0.03, 0.25])
cbar3 = plt.colorbar(p3, cax=cbax3, orientation='vertical', ticks=levs2[::2])
cbax3.tick_params(labelsize=fstic)
cbax3.set_ylabel("Fraction", fontsize=fslab)



#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig5.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig5_trans.png', dpi=300, bbox_inches='tight', transparent=True)


