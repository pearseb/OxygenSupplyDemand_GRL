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

data = nc.Dataset("ETOPO_JRA55_ndep_oxy200trend.nc",'r')
oxy_75s = data.variables['OXY200_75S'][...]
oxy_80s = data.variables['OXY200_80S'][...]

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()



#%% figure specifics

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.balance, 0.8)
levs1 = np.arange(-40, 41, 4)*0.1


contcol = 'black'
contwid=0.75

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)



#%% figure

fig = plt.figure(figsize=(8,8))
gs = GridSpec(2,1)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, oxy_75s, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
c1 = plt.contour(lons, lats, oxy_75s, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax2 = plt.subplot(gs[1,0], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, oxy_80s, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
c2 = plt.contour(lons, lats, oxy_80s, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

fig.subplots_adjust(top=0.825, bottom=0.025, left=0.05, right=0.95, hspace=0.15)


xx = 0.5; yy = 1.05
plt.text(xx, yy, '2005-2014 minus 1975-1984', fontsize=fslab, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, '2005-2014 minus 1980-1989', fontsize=fslab, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')

xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.2, 0.88, 0.6, 0.03])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='horizontal', ticks=levs1[::2])
cbax1.tick_params(labelsize=fstic, top=True, labeltop=True, bottom=False, labelbottom=False)
cbax1.set_xlabel("$\Delta$ O$_2$ ($\mu$M decade$^{-1}$)", fontsize=fslab)
cbax1.xaxis.set_label_position('top')


#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig-oxy200.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig-oxy200_trans.png', dpi=300, bbox_inches='tight', transparent=True)


