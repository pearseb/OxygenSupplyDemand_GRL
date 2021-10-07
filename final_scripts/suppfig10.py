# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:45:34 2021

Purpose
-------
    Figure of linear oxygen, TOU and ideal age trends averaged over 26.5-27.0 and
    27.0-27.4 potential density bounds 
    (hatching represents where the signal exceeds interannual noise)
    
    
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

# for colouring scatter plot according to density of data
from scipy.stats import gaussian_kde

# for lowess smoothing
from statsmodels.nonparametric.smoothers_lowess import lowess


#%% get data

os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//data_for_figures")


data = nc.Dataset("supfigure10.nc", 'r')
gmoc = data.variables['GMOC'][...]
amoc = data.variables['AMOC'][...]
pmoc = data.variables['GMOC'][...]
imoc = data.variables['GMOC'][...]
gwma = data.variables['WMA_GLO'][...]
awma = data.variables['WMA_ATL'][...]
pwma = data.variables['WMA_PAC'][...]
iwma = data.variables['WMA_IND'][...]

dep = data.variables['DEPTHT'][...]
lat = data.variables['ETOPO60Y'][...]

data.close()


#%% get differences in MOC and age

gmoc_ave = np.ma.average(gmoc[3:42,:,:],axis=0)
amoc_ave = np.ma.average(amoc[3:42,:,:],axis=0)
pmoc_ave = np.ma.average(pmoc[3:42,:,:],axis=0)
imoc_ave = np.ma.average(imoc[3:42,:,:],axis=0)

gwma_ave = np.ma.average(gwma[3:42,:,:],axis=0)
awma_ave = np.ma.average(awma[3:42,:,:],axis=0)
pwma_ave = np.ma.average(pwma[3:42,:,:],axis=0)
iwma_ave = np.ma.average(iwma[3:42,:,:],axis=0)

gmoc_dif = (np.ma.average(gmoc[33:42,:,:],axis=0) - np.ma.average(gmoc[3:12,:,:],axis=0))/30.*10
amoc_dif = (np.ma.average(amoc[33:42,:,:],axis=0) - np.ma.average(amoc[3:12,:,:],axis=0))/30.*10
pmoc_dif = (np.ma.average(pmoc[33:42,:,:],axis=0) - np.ma.average(pmoc[3:12,:,:],axis=0))/30.*10
imoc_dif = (np.ma.average(imoc[33:42,:,:],axis=0) - np.ma.average(imoc[3:12,:,:],axis=0))/30.*10

gwma_dif = (np.ma.average(gwma[33:42,:,:],axis=0) - np.ma.average(gwma[3:12,:,:],axis=0))/30.*10
awma_dif = (np.ma.average(awma[33:42,:,:],axis=0) - np.ma.average(awma[3:12,:,:],axis=0))/30.*10
pwma_dif = (np.ma.average(pwma[33:42,:,:],axis=0) - np.ma.average(pwma[3:12,:,:],axis=0))/30.*10
iwma_dif = (np.ma.average(iwma[33:42,:,:],axis=0) - np.ma.average(iwma[3:12,:,:],axis=0))/30.*10


#%% figure specifics 

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap5 = cmocean.tools.lighten(cmo.balance, 0.8)
levs5 = np.arange(-20,21,2)
cont5 = np.arange(0,1001,50)

colmap6 = cmocean.tools.lighten(cmo.balance, 0.8)
levs6 = np.arange(-10,11,2)*0.1
cont6 = np.array([-10,-5,-2,-1,-0.5,-0.2,0.0,0.2,0.5,1,2,5,10])


contcol = 'black'
contwid = 0.75
alf = 0.8


lat_labs1 = ['80$^{\circ}$S', '60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq ']


#%% make figure 


fig = plt.figure(figsize=(12,3.5))
gs = GridSpec(1,2)


lats, deps = np.meshgrid(lat,dep)

ax1 = plt.subplot(gs[0,0])
ax1.tick_params(labelsize=fstic)
p1 = plt.contourf(lats,deps, gmoc_ave, cmap=colmap5, levels=levs5, vmin=np.min(levs5), vmax=np.max(levs5), extend='both')
c1 = plt.contour(lats,deps, gwma_ave, colors=contcol, linewidths=contwid, levels=cont5)
plt.ylim(1200,10)
plt.xlim(-80,0)
plt.xticks(np.arange(-80,1,20), lat_labs1)

ax2 = plt.subplot(gs[0,1])
ax2.tick_params(labelsize=fstic, labelleft=False)
p2 = plt.contourf(lats,deps, gmoc_dif, cmap=colmap6, levels=levs6, vmin=np.min(levs6), vmax=np.max(levs6), extend='both')
c2 = plt.contour(lats,deps, gwma_dif, colors=contcol, linewidths=contwid, levels=cont6)
plt.ylim(1200,10)
plt.xlim(-80,0)
plt.xticks(np.arange(-80,1,20), lat_labs1)



fig.subplots_adjust(top=0.9, bottom=0.1, left=0.15, right=0.875, wspace=0.1, hspace=0.1)


xx = 0.5; yy = 1.05
plt.text(xx, yy, 'Global Meridional Overturning', fontsize=fslab, transform=ax1.transAxes, va='center', ha='center')
plt.text(xx, yy, 'trend (2005-2014 minus 1975-1984)', fontsize=fslab, transform=ax2.transAxes, va='center', ha='center')


xx = 0.025; yy = 1.05
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')



cbax5 = fig.add_axes([0.075, 0.15, 0.025, 0.7])
cbar5 = plt.colorbar(p1, cax=cbax5, orientation='vertical', ticks=levs5[::2])
cbax5.set_ylabel('overturning (Sv)', fontsize=fstic)
cbax5.tick_params(labelsize=fstic, right=False, labelright=False, left=True, labelleft=True)
cbax5.yaxis.set_label_position('left')

cbax6 = fig.add_axes([0.89, 0.15, 0.025, 0.7])
cbar6 = plt.colorbar(p2, cax=cbax6, orientation='vertical', ticks=levs6[::2])
cbax6.set_ylabel('$\Delta$ (Sv decade$^{-1}$)', fontsize=fstic)
cbax6.tick_params(labelsize=fstic)



#%% save figure

plt.clabel(c1, fmt='%i', manual=True, fontsize=fstic-2)
plt.clabel(c2, fmt='%.1f', manual=True, fontsize=fstic-2)


os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-suppfig10.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig10_trans.png', dpi=300, bbox_inches='tight', transparent=True)

