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



#%% get data

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")

data = nc.Dataset("figure_compare_ito_jra55_ndep.nc", 'r')
oxy_dep_dec = np.squeeze(data.variables['OXY_M_DEC'][...])
tou_dep_dec = np.squeeze(data.variables['TOU_M_DEC'][...])*(-1)

data = nc.Dataset("figure_compare_ito_jra55.nc", 'r')
oxy_pic_dec = np.squeeze(data.variables['OXY_M_DEC'][...])
tou_pic_dec = np.squeeze(data.variables['TOU_M_DEC'][...])*(-1)

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]

data = nc.Dataset('figure_ndep_important.nc', 'r')
ndep_important = data.variables['NDEP_IMPORTANT'][...]

data = nc.Dataset("figure_compare_ito_jra55_6thpanel.nc", 'r')
sol_trop = data.variables['SOL_TROP'][3:43]
sol_extra = data.variables['SOL_EXTRA'][3:43]
sol_glob = data.variables['SOL_GLOB'][3:43]
oxy_pic_trop = data.variables['OXY_PIC_TROP'][3:43]
oxy_pic_extra = data.variables['OXY_PIC_EXTRA'][3:43]
oxy_pic_glob = data.variables['OXY_PIC_GLOB'][3:43]
aou_pic_trop = data.variables['AOU_PIC_TROP'][3:43]*(-1)
aou_pic_extra = data.variables['AOU_PIC_EXTRA'][3:43]*(-1)
aou_pic_glob = data.variables['AOU_PIC_GLOB'][3:43]*(-1)
oxy_dep_trop = data.variables['OXY_DEP_TROP'][3:43]
oxy_dep_extra = data.variables['OXY_DEP_EXTRA'][3:43]
oxy_dep_glob = data.variables['OXY_DEP_GLOB'][3:43]
aou_dep_trop = data.variables['AOU_DEP_TROP'][3:43]*(-1)
aou_dep_extra = data.variables['AOU_DEP_EXTRA'][3:43]*(-1)
aou_dep_glob = data.variables['AOU_DEP_GLOB'][3:43]*(-1)

years = np.arange(1975,2015,1)

data.close()


#%% figure specifics 

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.balance, 0.8)
levs1 = np.array([-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
cont1 = levs1[::2]

contcol = 'black'
contwid = 0.5

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)


sol_cont = np.arange(-50,51,1)
aou_cont = np.arange(-50,51,1)
sol_cont, aou_cont = np.meshgrid(sol_cont, aou_cont)
oxy_cont = sol_cont+aou_cont


#%% make figure

fig = plt.figure(figsize=(6,9))
gs = GridSpec(15,1)

ax1 = plt.subplot(gs[0:6], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, oxy_dep_dec-oxy_pic_dec, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
plt.contourf(lons,lats, ndep_important, transform=ccrs.PlateCarree(), colors='none', levels=[0.9,1.1], hatches=['......'],zorder=2)
c1 = plt.contour(lons, lats, oxy_dep_dec-oxy_pic_dec, transform=ccrs.PlateCarree(), colors=contcol, linewidths=contwid, levels=cont1, zorder=2)

aou_dep_change = aou_dep_glob - aou_dep_glob[0]
aou_pic_change = aou_pic_glob - aou_pic_glob[0]
sol_change =  sol_glob - sol_glob[0]

ax2 = plt.subplot(gs[6:9,0])
ax2.tick_params(labelsize=fstic, labelbottom=False)
plt.plot(aou_pic_change, sol_change, color='k', label='Preindustrial N$_r$ deposition')
plt.plot(aou_dep_change, sol_change, color='forestgreen', label='Anthropogenic N$_r$ deposition')
c1 = plt.contour(aou_cont, sol_cont, oxy_cont, levels=np.arange(-100,101,0.2), colors=contcol, linewidths=contwid)
plt.ylim(-0.2,0.2)
plt.xlim(-1.5,0.05)
plt.ylabel('$\Delta$ O$_2^{sol}$ ($\mu$M)', fontsize=fslab)
plt.legend(frameon=False, loc='upper left', ncol=1, labelspacing=0.05, bbox_to_anchor=(0.2,1.45))

aou_dep_change = aou_dep_trop - aou_dep_trop[0]
aou_pic_change = aou_pic_trop - aou_pic_trop[0]
sol_change =  sol_trop - sol_trop[0]

ax3 = plt.subplot(gs[9:12,0])
ax3.tick_params(labelsize=fstic, labelbottom=False)
plt.plot(aou_pic_change, sol_change, color='k', label='Preindustrial N$_r$ deposition')
plt.plot(aou_dep_change, sol_change, color='forestgreen', label='Anthropogenic N$_r$ deposition')
c2 = plt.contour(aou_cont, sol_cont, oxy_cont, levels=np.arange(-100,101,0.2), colors=contcol, linewidths=contwid)
plt.ylim(-0.25,0.25)
plt.xlim(-1.5,0.05)
plt.ylabel('$\Delta$ O$_2^{sol}$ ($\mu$M)', fontsize=fslab)

aou_dep_change = aou_dep_extra - aou_dep_extra[0]
aou_pic_change = aou_pic_extra - aou_pic_extra[0]
sol_change =  sol_extra - sol_extra[0]

ax4 = plt.subplot(gs[12:15,0])
ax4.tick_params(labelsize=fstic, labelbottom=True)
plt.plot(aou_pic_change, sol_change, color='k', label='Preindustrial N$_r$ deposition')
plt.plot(aou_dep_change, sol_change, color='forestgreen', label='Anthropogenic N$_r$ deposition')
c3 = plt.contour(aou_cont, sol_cont, oxy_cont, levels=np.arange(-100,101,0.2), colors=contcol, linewidths=contwid)
plt.ylim(-0.25,0.25)
plt.xlim(-1.5,0.05)
plt.ylabel('$\Delta$ O$_2^{sol}$ ($\mu$M)', fontsize=fslab)
plt.xlabel('$\Delta$ O$_2^{AOU}$ ($\mu$M)', fontsize=fslab)


fig.subplots_adjust(top=0.875, bottom=0.1, left=0.15, right=0.9, hspace=1.5)


xx = 0.5; yy = 1.05
plt.text(xx, yy, '$\Delta$ O$_2$ due to N$_r$ deposition', fontsize=fstic, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')

xx = 1.05; yy = 0.5
plt.text(xx,yy,'Global', transform=ax2.transAxes, rotation=90, ha='center', va='center', fontsize=fstic)
plt.text(xx,yy,'20$^{\circ}$S-20$^{\circ}$N', transform=ax3.transAxes, rotation=90, ha='center', va='center', fontsize=fstic)
plt.text(xx,yy,'Extra-tropics', transform=ax4.transAxes, rotation=90, ha='center', va='center', fontsize=fstic)

xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx-0.025, yy+0.15, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx-0.025, yy+0.15, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx-0.025, yy+0.15, 'd', fontsize=fslab+2, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.2, 0.9, 0.65, 0.025])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='horizontal', ticks=levs1[::2])
cbax1.set_xlabel('$\Delta$ O$_2$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
cbax1.tick_params(labelsize=fstic, top=True, labeltop=True, bottom=False, labelbottom=False)
cbax1.xaxis.set_label_position('top')


#%% save figure

plt.clabel(c1, manual=True, fmt='%.1f', fontsize=fstic-2)
plt.clabel(c2, manual=True, fmt='%.1f', fontsize=fstic-2)
plt.clabel(c3, manual=True, fmt='%.1f', fontsize=fstic-2)


os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig1.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig1_trans.png', dpi=300, bbox_inches='tight', transparent=True)

