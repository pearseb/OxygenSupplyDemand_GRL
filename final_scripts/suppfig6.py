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
oxy_tre = data.variables['OXY_TRE_ZAVE'][...]
tou_tre = data.variables['TOU_TRE_ZAVE'][...]
wma_tre = data.variables['WMA_TRE_ZAVE'][...]
rem_tre = data.variables['RDEM_TRE_ZAVE'][...]
wma_r = data.variables['WMA_R'][...]
rem_r = data.variables['RDEM_R'][...]
wma_s = data.variables['WMA_SLOPE'][...]
rem_s = data.variables['RDEM_SLOPE'][...]
wma_std = data.variables['WMA_STD'][...]
rem_std = data.variables['RDEM_STD'][...]

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()



#%% figure specifics

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap = cmocean.tools.lighten(cmo.balance, 0.8)
cont1 = np.array([-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5])
cont2 = np.arange(-5,5.1,0.5)
cont3 = cont1 * 0.1


contcol = 'black'
contwid=0.75

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)



#%% figure

fig = plt.figure(figsize=(14,9.5))
gs = GridSpec(3,2)

ax1 = plt.subplot(gs[0,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, oxy_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=cont1, vmin=np.min(cont1),vmax=np.max(cont1), zorder=1, extend='both')
c1 = plt.contour(lons, lats, oxy_tre, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax2 = plt.subplot(gs[0,1], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, tou_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=cont1, vmin=np.min(cont1),vmax=np.max(cont1), zorder=1, extend='both')
c2 = plt.contour(lons, lats, tou_tre, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax3 = plt.subplot(gs[1,0], projection=proj)
ax3.tick_params(labelsize=fstic)
ax3.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax3.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax3.coastlines()
p3 = plt.contourf(lons, lats, wma_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=cont2, vmin=np.min(cont2),vmax=np.max(cont2), zorder=1, extend='both')
c3 = plt.contour(lons, lats, wma_tre, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)

ax4 = plt.subplot(gs[1,1], projection=proj)
ax4.tick_params(labelsize=fstic)
ax4.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax4.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax4.coastlines()
p4 = plt.contourf(lons, lats, rem_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=cont3, vmin=np.min(cont3),vmax=np.max(cont3), zorder=1, extend='both')
c4 = plt.contour(lons, lats, rem_tre, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)


ax5 = plt.subplot(gs[2,0], projection=proj)
ax5.tick_params(labelsize=fstic)
ax5.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax5.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax5.coastlines()
p5 = plt.contourf(lons, lats, wma_s*wma_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=cont1, vmin=np.min(cont1),vmax=np.max(cont1), zorder=1, extend='both')
c5 = plt.contour(lons, lats, wma_s*wma_tre, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)
plt.contourf(lons, lats, wma_r, transform=ccrs.PlateCarree(), colors='none', levels=[-1,0,1], hatches=['......',' '], zorder=2)

ax6 = plt.subplot(gs[2,1], projection=proj)
ax6.tick_params(labelsize=fstic)
ax6.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax6.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax6.coastlines()
p6 = plt.contourf(lons, lats, rem_s*rem_tre, transform=ccrs.PlateCarree(), cmap=colmap, levels=cont1, vmin=np.min(cont1),vmax=np.max(cont1), zorder=1, extend='both')
c6 = plt.contour(lons, lats, rem_s*rem_tre, transform=ccrs.PlateCarree(), levels=[0], zorder=2, colors=contcol, linewidths=contwid)
plt.contourf(lons, lats, rem_r, transform=ccrs.PlateCarree(), colors='none', levels=[-1,0,1], hatches=['......',' '], zorder=2)

fig.subplots_adjust(top=0.95, bottom=0.05, left=0.15, right=0.85, wspace=0.05, hspace=0.05)


xx = 0.5; yy = 1.05
plt.text(xx, yy, '$\Delta$ O$_2$', fontsize=fslab, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, '$\Delta$ AOU', fontsize=fslab, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, '$\Delta$ Age', fontsize=fslab, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, '$\Delta$ demand', fontsize=fslab, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy+0.025, r'$\frac{\Delta AOU}{\Delta Age}$', fontsize=fslab, transform=ax5.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy+0.025, r'$\frac{\Delta AOU}{\Delta demand}$', fontsize=fslab, transform=ax6.transAxes, va='center', ha='center', fontweight='bold')

xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'd', fontsize=fslab+2, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'e', fontsize=fslab+2, transform=ax5.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'f', fontsize=fslab+2, transform=ax6.transAxes, va='center', ha='center', fontweight='bold')



cbax1 = fig.add_axes([0.09, 0.68, 0.03, 0.25])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=cont1[::3])
cbax1.tick_params(labelsize=fstic, left=True, labelleft=True, right=False, labelright=False)
cbax1.set_ylabel("$\mu$M decade$^{-1}$", fontsize=fslab)
cbax1.yaxis.set_label_position('left')

cbax3 = fig.add_axes([0.09, 0.375, 0.03, 0.25])
cbar3 = plt.colorbar(p3, cax=cbax3, orientation='vertical', ticks=cont2[::2])
cbax3.tick_params(labelsize=fstic, left=True, labelleft=True, right=False, labelright=False)
cbax3.set_ylabel("years decade$^{-1}$", fontsize=fslab)
cbax3.yaxis.set_label_position('left')

cbax5 = fig.add_axes([0.09, 0.07, 0.03, 0.25])
cbar5 = plt.colorbar(p5, cax=cbax5, orientation='vertical', ticks=cont1[::3])
cbax5.tick_params(labelsize=fstic, left=True, labelleft=True, right=False, labelright=False)
cbax5.set_ylabel("$\mu$M decade$^{-1}$", fontsize=fslab)
cbax5.yaxis.set_label_position('left')


cbax2 = fig.add_axes([0.88, 0.68, 0.03, 0.25])
cbar2 = plt.colorbar(p2, cax=cbax2, orientation='vertical', ticks=cont1[::3])
cbax2.tick_params(labelsize=fstic)
cbax2.set_ylabel('$\mu$M decade$^{-1}$', fontsize=fslab)

cbax4 = fig.add_axes([0.88, 0.375, 0.03, 0.25])
cbar4 = plt.colorbar(p4, cax=cbax4, orientation='vertical', ticks=cont3[::3])
cbax4.tick_params(labelsize=fstic)
cbax4.set_ylabel('$\mu$M decade$^{-1}$', fontsize=fslab)

cbax6 = fig.add_axes([0.88, 0.07, 0.03, 0.25])
cbar6 = plt.colorbar(p6, cax=cbax6, orientation='vertical', ticks=cont1[::3])
cbax6.tick_params(labelsize=fstic)
cbax6.set_ylabel('$\mu$M decade$^{-1}$', fontsize=fslab)



#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig6.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig6_trans.png', dpi=300, bbox_inches='tight', transparent=True)


