# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 10:39:44 2021

Purpose
-------
    Create figure of change in SST between 2014 and 1984
    

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

data = nc.Dataset("suppfigure_sst_change.nc",'r')
obs = np.ma.squeeze(data.variables['OBS_SST'][...])
mod = np.ma.squeeze(data.variables['MOD_SST'][...])

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()


#%% figure specifics

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap = cmocean.tools.lighten(cmo.balance, 0.8)
levs = np.array([-2, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
cont = [0]
contcol = 'black'
contwid=0.75

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)


#%% figure

fig = plt.figure(figsize=(8,9))
gs = GridSpec(2,1)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, obs, transform=ccrs.PlateCarree(), cmap=colmap, levels=levs, vmin=np.min(levs),vmax=np.max(levs), zorder=1, extend='both')
#c1 = plt.contour(lons, lats, obs, transform=ccrs.PlateCarree(), levels=cont, zorder=2, colors=contcol, linewidths=contwid)

ax2 = plt.subplot(gs[1,0], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, mod, transform=ccrs.PlateCarree(), cmap=colmap, levels=levs, vmin=np.min(levs),vmax=np.max(levs), zorder=1, extend='both')
#c2 = plt.contour(lons, lats, mod, transform=ccrs.PlateCarree(), levels=cont, zorder=2, colors=contcol, linewidths=contwid)


fig.subplots_adjust(top=0.83, bottom=0.03, left=0.05, right=0.95, wspace=0.05, hspace=0.15)


xx = 0.5; yy = 1.05
plt.text(xx, yy, 'NOAA OI-SST v2', fontsize=fslab, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'Hindcast simulation', fontsize=fslab, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')

xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')

cbax1 = fig.add_axes([0.2, 0.885, 0.6, 0.03])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='horizontal', ticks=levs[::2])
cbax1.tick_params(labelsize=fstic, bottom=False, labelbottom=False, top=True, labeltop=True)
cbax1.set_xlabel("$\Delta$ SST ($^{\circ}$C) 2014 minus 1984", fontsize=fslab, labelpad=10)
cbax1.xaxis.set_label_position('top')


#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig8.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig8_trans.png', dpi=300, bbox_inches='tight', transparent=True)



