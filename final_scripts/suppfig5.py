# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:45:34 2021

Purpose
-------
    Figure of linear oxygen, solubility and AOU trends averaged between 200-1000 metres
    in the hindcast model simulation and in the Helm climatology from 1970-1992
    
    
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

from scipy.optimize import curve_fit


#%% get data

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")

data = nc.Dataset("figure3.nc", 'r')
ver = np.ma.squeeze(data.variables['WO_TRE'][...])      # vertical velocity (cm/day per decade)

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()


#%% figure specifics 

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.balance, 0.8)
levs1 = np.arange(-20,20.1,2)*0.1
cont1 = levs1[::5]


contcol = 'black'
contwid = 0.5


proj = ccrs.Geostationary(central_longitude=-15, satellite_height=2.5e7)
#proj = ccrs.Orthographic(central_longitude=-40, central_latitude=0.0)
lons,lats = np.meshgrid(lon,lat)


#%% make figure

fig = plt.figure(figsize=(8,6))
gs = GridSpec(1,1)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
#p1 = plt.pcolormesh(lons, lats, ver, transform=ccrs.PlateCarree(), cmap=colmap1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1)
p1 = plt.contourf(lons, lats, ver, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')


fig.subplots_adjust(top=0.925, bottom=0.05, left=0.055, right=0.80)


xx = 0.5; yy = 1.05
plt.text(xx, yy, 'Change in upwelling intensity at 50 metres', fontsize=fstic, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.81, 0.2, 0.035, 0.6])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=levs1[::2])
cbax1.set_ylabel('$\Delta$ cm day$^{-1}$ per decade', fontsize=fslab)
cbax1.tick_params(labelsize=fstic)


#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig5.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig5.png', dpi=300, bbox_inches='tight', transparent=True)

