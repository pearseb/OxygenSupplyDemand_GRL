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
npp = np.ma.squeeze(data.variables['NPP_TRE'][...])     # net primary production (gC/m2 per decade)
bex = np.ma.squeeze(data.variables['BEX_TRE'][...])     # b exponent of Martin Curve (per decade)

atlnpp_rbio = np.ma.squeeze(data.variables['NPP_TRE_ATLMASK'][...])     # net primary production (gC/m2 per decade)
atlbex_rbio = np.ma.squeeze(data.variables['BEX_TRE_ATLMASK'][...])     # b exponent of Martin Curve (per decade)
atldem_rbio = np.ma.squeeze(data.variables['RDEM_TRE_ATLMASK'][...])*800    # biological oxygen demand (mmol O2 /m2 per decade)
atlmask = np.ma.squeeze(data.variables['ATLMASK'][...])    # mask of grid cells where scatter plot points were taken from

pacnpp_rbio = np.ma.squeeze(data.variables['NPP_TRE_PACMASK'][...])     # net primary production (gC/m2 per decade)
pacbex_rbio = np.ma.squeeze(data.variables['BEX_TRE_PACMASK'][...])     # b exponent of Martin Curve (per decade)
pacdem_rbio = np.ma.squeeze(data.variables['RDEM_TRE_PACMASK'][...])*800    # biological oxygen demand (mmol O2 /m2 per decade)
pacmask = np.ma.squeeze(data.variables['PACMASK'][...])    # mask of grid cells where scatter plot points were taken from


lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]


data.close()


#%% figure specifics 

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.delta, 0.8)
levs1 = np.arange(-30,30.1,3)
cont1 = levs1[::5]

colmap2 = cmocean.tools.lighten(cmo.balance, 0.8)
levs2 = np.arange(-5,5.1,0.5)*0.01
cont2 = levs2[::5]

contcol = 'black'
contwid = 0.5


alf = 0.2
col = 'grey'
mar = 'o'
edc = 'k'
liw = 0.5
zor = 1

x1a=-50;x2a=50
y1a=-0.3*800;y2a=0.3*800

x1b=-0.15;x2b=0.15
y1b=-0.3*800;y2b=0.3*800


zerocol='k'
zerosty='--'
zerowid=0.75
zerozor=5

onecol='k'
onesty='-'
onewid=1.5
onezor=5

proj = ccrs.Geostationary(central_longitude=-60, satellite_height=3.5e7)
proj = ccrs.Orthographic(central_longitude=-60, central_latitude=0.0)
proj = ccrs.Robinson(central_longitude=-60)
lons,lats = np.meshgrid(lon,lat)


#%% make figure

fig = plt.figure(figsize=(10,9.5))
gs = GridSpec(15,2)

ax1 = plt.subplot(gs[0:6,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=3, facecolor='silver')
ax1.coastlines(zorder=3)
ax1.set_extent([-140,40,-50,50], crs=ccrs.PlateCarree())
plt.contourf(lons, lats, atlmask, transform=ccrs.PlateCarree(), colors='none', hatches=['....'], levels=[0.5,1.5], zorder=1)
plt.contourf(lons, lats, pacmask, transform=ccrs.PlateCarree(), colors='none', hatches=['xxxx'], levels=[0.5,1.5], zorder=1)
p1 = plt.contourf(lons, lats, npp, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
c1 = plt.contour(lons, lats, npp, transform=ccrs.PlateCarree(), colors=contcol, linewidths=contwid, levels=cont1, zorder=2)


#plt.text(-12,15, "Canary", transform=ccrs.PlateCarree(), fontsize=fslab, va='center', ha='left')
#plt.text(15,-10, "Benguela", transform=ccrs.PlateCarree(), fontsize=fslab, va='center', ha='left', rotation=75)


ax2 = plt.subplot(gs[0:6,1], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=3, facecolor='silver')
ax2.coastlines(zorder=3)
ax2.set_extent([-140,40,-50,50], crs=ccrs.PlateCarree())
p2 = plt.contourf(lons, lats, bex, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='both')
plt.contour(lons, lats, bex, transform=ccrs.PlateCarree(), colors=contcol, linewidths=contwid, levels=cont2, zorder=2)


ax3 = plt.subplot(gs[6:10,0])
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.tick_params(labelsize=fstic)
#plt.scatter(np.ma.compressed(npp_rbio), np.ma.compressed(dem_rbio), color=col, alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)
plt.scatter(np.ma.compressed(atlnpp_rbio), np.ma.compressed(atldem_rbio), c=np.ma.compressed(atlbex_rbio), \
            cmap=colmap2, vmin=np.min(levs2), vmax=np.max(levs2), alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)


def linear(x,a,b):
    return a*x + b

popt,pcov = curve_fit(linear, np.ma.compressed(atlnpp_rbio), np.ma.compressed(atldem_rbio))
r = np.corrcoef(np.ma.compressed(atlnpp_rbio),np.ma.compressed(atldem_rbio))[0,1]
plt.plot(np.linspace(x1a,x2a,100),linear(np.linspace(x1a,x2a,100), popt[0], popt[1]), color='firebrick', linestyle=onesty, linewidth=onewid, zorder=onezor)

plt.text(0.25,0.85, 'r = %.2f'%(r), fontsize=fslab, transform=ax3.transAxes, va='center', ha='center')
plt.text(0.25,0.75, 'slope = %.2f'%(popt[0]), fontsize=fslab, transform=ax3.transAxes, va='center', ha='center')

plt.plot((x1a,x2a),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1a,y2a), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1a,y2a); plt.xlim(x1a,x2a)
plt.ylabel('$\Delta$ demand per decade\n(mmol O$_2$ m$^{-2}$)', fontsize=fslab)


ax4 = plt.subplot(gs[6:10,1])
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.tick_params(labelsize=fstic, labelleft=False)
#plt.scatter(np.ma.compressed(atlbex_rbio), np.ma.compressed(atldem_rbio), color=col, alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)
plt.scatter(np.ma.compressed(atlbex_rbio), np.ma.compressed(atldem_rbio), c=np.ma.compressed(atlnpp_rbio), \
            cmap=colmap1, vmin=np.min(levs1), vmax=np.max(levs1), alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)

def linear(x,a,b):
    return a*x + b

popt,pcov = curve_fit(linear, np.ma.compressed(atlbex_rbio), np.ma.compressed(atldem_rbio))
r = np.corrcoef(np.ma.compressed(atlbex_rbio),np.ma.compressed(atldem_rbio))[0,1]
plt.plot(np.linspace(x1a,x2a,100),linear(np.linspace(x1a,x2a,100), popt[0], popt[1]), color='firebrick', linestyle=onesty, linewidth=onewid, zorder=onezor)

plt.text(0.25,0.85, 'r = %.2f'%(r), fontsize=fslab, transform=ax4.transAxes, va='center', ha='center')
plt.text(0.25,0.75, 'slope = %i'%(popt[0]), fontsize=fslab, transform=ax4.transAxes, va='center', ha='center')

plt.plot((x1b,x2b),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1b,y2b), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1b,y2b); plt.xlim(x1b,x2b)
plt.xticks(np.arange(-10,15.1,5)*0.01, np.arange(-10,15.1,5)*0.01)


ax5 = plt.subplot(gs[11:15,0])
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.tick_params(labelsize=fstic)
plt.scatter(np.ma.compressed(pacnpp_rbio), np.ma.compressed(pacdem_rbio), c=np.ma.compressed(pacbex_rbio), \
            cmap=colmap2, vmin=np.min(levs2), vmax=np.max(levs2), alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)


def linear(x,a,b):
    return a*x + b

popt,pcov = curve_fit(linear, np.ma.compressed(pacnpp_rbio), np.ma.compressed(pacdem_rbio))
r = np.corrcoef(np.ma.compressed(pacnpp_rbio),np.ma.compressed(pacdem_rbio))[0,1]
plt.plot(np.linspace(x1a,x2a,100),linear(np.linspace(x1a,x2a,100), popt[0], popt[1]), color='firebrick', linestyle=onesty, linewidth=onewid, zorder=onezor)

plt.text(0.25,0.85, 'r = %.2f'%(r), fontsize=fslab, transform=ax5.transAxes, va='center', ha='center')
plt.text(0.25,0.75, 'slope = %.2f'%(popt[0]), fontsize=fslab, transform=ax5.transAxes, va='center', ha='center')

plt.plot((x1a,x2a),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1a,y2a), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1a,y2a); plt.xlim(x1a,x2a)
plt.ylabel('$\Delta$ demand per decade\n(mmol O$_2$ m$^{-2}$)', fontsize=fslab)
plt.xlabel('$\Delta$ NPP per decade (g C m$^{-2}$)', fontsize=fslab)


ax6 = plt.subplot(gs[11:15,1])
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.tick_params(labelsize=fstic, labelleft=False)
plt.scatter(np.ma.compressed(pacbex_rbio), np.ma.compressed(pacdem_rbio), c=np.ma.compressed(pacnpp_rbio), \
            cmap=colmap1, vmin=np.min(levs1), vmax=np.max(levs1), alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)

def linear(x,a,b):
    return a*x + b

popt,pcov = curve_fit(linear, np.ma.compressed(pacbex_rbio), np.ma.compressed(pacdem_rbio))
r = np.corrcoef(np.ma.compressed(pacbex_rbio),np.ma.compressed(pacdem_rbio))[0,1]
plt.plot(np.linspace(x1a,x2a,100),linear(np.linspace(x1a,x2a,100), popt[0], popt[1]), color='firebrick', linestyle=onesty, linewidth=onewid, zorder=onezor)

plt.text(0.25,0.85, 'r = %.2f'%(r), fontsize=fslab, transform=ax6.transAxes, va='center', ha='center')
plt.text(0.25,0.75, 'slope = %i'%(popt[0]), fontsize=fslab, transform=ax6.transAxes, va='center', ha='center')

plt.plot((x1b,x2b),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1b,y2b), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1b,y2b); plt.xlim(x1b,x2b)
plt.xlabel("$\Delta$ Martin's exponent per decade", fontsize=fslab)
plt.xticks(np.arange(-10,15.1,5)*0.01, np.arange(-10,15.1,5)*0.01)



fig.subplots_adjust(top=0.925, bottom=0.075, left=0.125, right=0.95, wspace=0.05)



xx = 0.05; yy = 1.05
plt.text(xx, yy+0.2, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy+0.2, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'd', fontsize=fslab+2, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'e', fontsize=fslab+2, transform=ax5.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'f', fontsize=fslab+2, transform=ax6.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.175, 0.9, 0.3, 0.03])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='horizontal', ticks=levs1[::4])
cbax1.set_xlabel('$\Delta$ NPP per decade (g C m$^{-2}$)', fontsize=fslab)
cbax1.tick_params(labelsize=fstic, top=True, labeltop=True, bottom=False, labelbottom=False)
cbax1.xaxis.set_label_position('top')

cbax2 = fig.add_axes([0.6, 0.9, 0.3, 0.03])
cbar2 = plt.colorbar(p2, cax=cbax2, orientation='horizontal', ticks=levs2[::4])
cbax2.set_xlabel("$\Delta$ Martin's exponent per decade", fontsize=fslab)
cbax2.tick_params(labelsize=fstic, top=True, labeltop=True, bottom=False, labelbottom=False)
cbax2.xaxis.set_label_position('top')

plt.text(1.12,0.5,'Deeper', transform=cbax2.transAxes, ha='center', va='center', fontsize=fstic)
plt.text(-0.15,0.5,'Shallower', transform=cbax2.transAxes, ha='center', va='center', fontsize=fstic)

xx = 0.5; yy = 1.05
plt.text(xx,yy,'Atlantic', transform=ax3.transAxes, ha='center', va='center', fontsize=fstic, fontweight='bold')
plt.text(xx,yy,'Atlantic', transform=ax4.transAxes, ha='center', va='center', fontsize=fstic, fontweight='bold')
plt.text(xx,yy,'Pacific', transform=ax5.transAxes, ha='center', va='center', fontsize=fstic, fontweight='bold')
plt.text(xx,yy,'Pacific', transform=ax6.transAxes, ha='center', va='center', fontsize=fstic, fontweight='bold')


#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-fig3.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-fig3.eps', dpi=300, bbox_inches='tight')
fig.savefig('fig-fig3_trans.png', dpi=300, bbox_inches='tight', transparent=True)

