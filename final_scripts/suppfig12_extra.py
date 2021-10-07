# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:45:34 2021

Purpose
-------
    
    
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

tou_jra55 = np.ma.squeeze(data.variables['TOU_TRE_MASK'][...])*800     # TOU (mmol O2 / m2 per decade)
npp_jra55 = np.ma.squeeze(data.variables['NPP_TRE_MASK'][...])     # net primary production (gC/m2 per decade)
bex_jra55 = np.ma.squeeze(data.variables['BEX_TRE_MASK'][...])     # b exponent of Martin Curve (per decade)
dem_jra55 = np.ma.squeeze(data.variables['RDEM_TRE_MASK'][...])*800    # biological oxygen demand (mmol O2 /m2 per decade)

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]

data.close()


#%% figure specifics 

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

contcol = 'black'
contwid = 0.5


alf = 0.2
col = 'grey'
mar = 'o'
edc = 'k'
liw = 0.5
zor = 1

cols=['royalblue', 'goldenrod', 'firebrick']


x1a=-50;x2a=50
y1a=-250;y2a=250

x1b=-250;x2b=250
y1b=-1500;y2b=1500


zerocol='k'
zerosty='--'
zerowid=0.75
zerozor=5

onecol='k'
onesty='-'
onewid=1.5
onezor=5


#%% make figure

fig = plt.figure(figsize=(8,8))
gs = GridSpec(1,1)

ax1 = plt.subplot(gs[0])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(labelsize=fstic)
plt.scatter(np.ma.compressed(npp_jra55), np.ma.compressed(dem_jra55), color=col, alpha=alf, zorder=zor, marker=mar, edgecolor=edc, linewidth=liw)

def linear(x,a,b):
    return a*x + b

popt,pcov = curve_fit(linear, np.ma.compressed(npp_jra55), np.ma.compressed(dem_jra55))
r2 = np.corrcoef(np.ma.compressed(npp_jra55),np.ma.compressed(dem_jra55))[0,1]**2
plt.plot(np.linspace(x1a,x2a,100),linear(np.linspace(x1a,x2a,100), popt[0], popt[1]), color='k', linestyle=onesty, linewidth=onewid, zorder=onezor)

plt.text(0.25,0.85, 'r$^2$ = %.2f'%(r2), fontsize=fslab, transform=ax1.transAxes, va='center', ha='center')
plt.text(0.25,0.75, 'slope = %.2f $\pm$ %.2f'%(popt[0], np.sqrt(np.diag(pcov))[0]), fontsize=fslab, transform=ax1.transAxes, va='center', ha='center')

plt.plot((x1a,x2a),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1a,y2a), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1a,y2a); plt.xlim(x1a,x2a)
plt.ylabel('$\Delta$ mesopelagic demand\n($\mu$M O$_2$ m$^{-2}$ decade$^{-1}$)', fontsize=fslab)
plt.xlabel('$\Delta$ depth-integrated NPP\n(g C m$^{-2}$ decade$^{-1}$)', fontsize=fslab)

plt.subplots_adjust(left=0.175, bottom=0.175)

#%% save figure

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig12_extra.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig12_extra_trans.png', dpi=300, bbox_inches='tight', transparent=True)
