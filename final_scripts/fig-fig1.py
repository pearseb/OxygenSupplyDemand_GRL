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
oxy_o_dec = np.squeeze(data.variables['OXY_O_DEC'][...])
oxy_m_dec = np.squeeze(data.variables['OXY_M_DEC'][...])
sol_o_dec = np.squeeze(data.variables['SOL_O_DEC'][...])
sol_m_dec = np.squeeze(data.variables['SOL_M_DEC'][...])
aou_o_dec = np.squeeze(data.variables['AOU_O_DEC'][...])*(-1)
aou_m_dec = np.squeeze(data.variables['AOU_M_DEC'][...])*(-1)
tou_m_dec = np.squeeze(data.variables['TOU_M_DEC'][...])*(-1)
dis_m_dec = np.squeeze(data.variables['DIS_M_DEC'][...])*(-1)
sup_m_dec = np.squeeze(data.variables['SUP_M_DEC'][...])*(-1)
dem_m_dec = np.squeeze(data.variables['DEM_M_DEC'][...])*(-1)

lon = data.variables['ETOPO60X'][...]
lat = data.variables['ETOPO60Y'][...]

'''
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
'''
years = np.arange(1975,2015,1)


data = nc.Dataset("cmip6_multimodel_trends.nc",'r')
obs = data.variables['OXY_ITO'][...]
access = data.variables['OXY_ACCESS'][...]
canesm = data.variables['OXY_CANESM'][...]
ipsl = data.variables['OXY_IPSL'][...]
miroc = data.variables['OXY_MIROC'][...]
mpi = data.variables['OXY_MPI'][...]
mri = data.variables['OXY_MRI'][...]
ukesm = data.variables['OXY_UKESM'][...]
jra = data.variables['OXY_JRA55'][0:62]

m_access = data.variables['OXY_OBMASK_ACCESS'][0:64]
m_canesm = data.variables['OXY_OBMASK_CANESM'][0:64]
m_ipsl = data.variables['OXY_OBMASK_IPSL'][0:64]
m_miroc = data.variables['OXY_OBMASK_MIROC'][0:64]
m_mpi = data.variables['OXY_OBMASK_MPI'][0:64]
m_mri = data.variables['OXY_OBMASK_MRI'][0:64]
m_ukesm = data.variables['OXY_OBMASK_UKESM'][0:64]


years_cmip6 = np.arange(1951,2016,1)
years_jra = np.arange(1958,2020,1)

data.close()


#%% get the range of CMIP6 and the linear trends of obs, the hindcast and CMIP6 multi-model mean

cmip6 = np.transpose(np.vstack((access, canesm, ipsl, miroc, mpi, mri, ukesm)))
m_cmip6 = np.transpose(np.vstack((m_access, m_canesm, m_ipsl, m_miroc, m_mpi, m_mri, m_ukesm)))


from scipy.optimize import curve_fit

def linear(x,a,b):
    return a*x + b

ito_popt,ito_pcov = curve_fit(linear, years_cmip6[24:64], obs[24:64])
cmip6_popt,cmip6_pcov = curve_fit(linear, years_cmip6[24:64], np.mean(cmip6[24:64],axis=1))
m_cmip6_popt,m_cmip6_pcov = curve_fit(linear, years_cmip6[24:64], np.mean(m_cmip6[24:64],axis=1))
jra_popt,jra_pcov = curve_fit(linear, years_jra[17:57], jra[17:57])

ito_error = np.sqrt(np.diag(ito_pcov))[0]
cmip6_error = np.sqrt(np.diag(cmip6_pcov))[0]
m_cmip6_error = np.sqrt(np.diag(m_cmip6_pcov))[0]
jra_error = np.sqrt(np.diag(jra_pcov))[0]

print(ito_popt[0], ito_error)
print(cmip6_popt[0], cmip6_error)
print(jra_popt[0], jra_error)


#%% figure specifics 

fslab = 15
fstic = 13
lw = 0.75
gridalf = 0.5

colmap1 = cmocean.tools.lighten(cmo.balance, 0.8)
levs1 = np.array([-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5])
cont1 = levs1[::5]

colmap2 = cmocean.tools.lighten(cmo.balance, 0.8)
levs2 = np.arange(-60,61,8)*0.1
cont2 = levs2[::5]


contcol = 'k'; contwid = 0.5; contsty = '--'; contalf = 0.5; contzor = 2

alf = [0.1, 0.1, 0.1, 0.1]
col = ['grey', 'firebrick', 'goldenrod', 'royalblue']
mar = ['o', 'o', 'o', 'o']
edc = ['k', 'k', 'k', 'k']
liw = [0.5, 0.5, 0.5, 0.5]
zor = [1, 2, 3, 4]

x1a=-10;x2a=10
y1a=-5;y2a=5

x1b=-20;x2b=20
y1b=-10;y2b=10

x1c=-10;x2c=10
y1c=-5;y2c=5


zerocol='k'
zerosty='--'
zerowid=0.75
zerozor=5

onecol='k'
onesty='-'
onewid=1.5
onezor=5

lab = ['O$_2^{AOU}$', 'O$_2^{sol}$', 'O$_2^{dis}$', 'O$_2^{TOU}$', 'O$_2^{sup}$', 'O$_2^{dem}$']

proj = ccrs.Robinson(central_longitude=200)
lons,lats = np.meshgrid(lon,lat)


sol_cont = np.arange(-50,51,1)
aou_cont = np.arange(-50,51,1)
sol_cont, aou_cont = np.meshgrid(sol_cont, aou_cont)
oxy_cont = sol_cont+aou_cont


### For panel f
fcol = ['k', 'firebrick', 'grey']
fwid = [2, 2, 2]
fsty = ['-', '-', '-']
falf = [1, 1, 0.75]
fzor = [1, 1, 1]
flab = ["Ito et al. (2017)    %.2f $\pm$ %.2f $\mu$M decade$^{-1}$ since 1975"%(ito_popt[0]*10, ito_error*10), \
       "hindcast               %.2f $\pm$ %.2f $\mu$M decade$^{-1}$ since 1975"%(jra_popt[0]*10, jra_error*10), \
       "CMIP6 range       %.2f $\pm$ %.2f $\mu$M decade$^{-1}$ since 1975 (multi-model mean)"%(cmip6_popt[0]*10, cmip6_error*10)]

flab = ["Ito et al. (2017)", "Hindcast simulation", "CMIP6 historical"]



#%% make figure

fig = plt.figure(figsize=(11,8))
gs = GridSpec(32,2)

ax1 = plt.subplot(gs[0:12,0], projection=proj)
ax1.tick_params(labelsize=fstic)
ax1.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax1.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax1.coastlines()
p1 = plt.contourf(lons, lats, oxy_m_dec, transform=ccrs.PlateCarree(), cmap=colmap1, levels=levs1, vmin=np.min(levs1),vmax=np.max(levs1), zorder=1, extend='both')
c1 = plt.contour(lons, lats, oxy_m_dec, transform=ccrs.PlateCarree(), colors=contcol, linewidths=contwid, levels=cont1, zorder=2)

ax2 = plt.subplot(gs[0:9,1])
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.tick_params(labelsize=fstic, left=False, labelleft=False, right=True, labelright=True)
plt.plot(years_cmip6, obs*0, color=contcol, linewidth=contwid, linestyle=contsty, alpha=contalf, zorder=contzor)
plt.plot(years_cmip6, obs, color=fcol[0], linewidth=fwid[0], linestyle=fsty[0], alpha=falf[0], zorder=fzor[0], label=flab[0])
plt.plot(years_jra, jra, color=fcol[1], linewidth=fwid[1], linestyle=fsty[1], alpha=falf[1], zorder=fzor[1], label=flab[1])
plt.plot(years_cmip6[0:64], np.min(cmip6,axis=1), color=fcol[2], linewidth=fwid[2], linestyle=fsty[2], alpha=falf[2], zorder=fzor[2])
plt.plot(years_cmip6[0:64], np.max(cmip6,axis=1), color=fcol[2], linewidth=fwid[2], linestyle=fsty[2], alpha=falf[2], zorder=fzor[2])
plt.fill_between(years_cmip6[0:64], np.min(cmip6,axis=1), np.max(cmip6,axis=1), color=fcol[2], alpha=falf[2], label=flab[2])
plt.legend(frameon=False, loc='lower left', ncol=1, fontsize=11, labelspacing=0.25, columnspacing=0.5)
plt.xlim(1951,2015)
plt.ylim(-6,3)
#plt.xlabel('Year', fontsize=fslab)
plt.ylabel('$\Delta$ O$_2$ ($\mu$M)', fontsize=fslab)
ax2.yaxis.set_label_position('right')

'''
ax2 = plt.subplot(gs[0:6,1], projection=proj)
ax2.tick_params(labelsize=fstic)
ax2.gridlines(linestyle='--', linewidth=0.5, color='grey', alpha=gridalf, zorder=3)
ax2.add_feature(cartopy.feature.LAND, zorder=1, facecolor='silver')
ax2.coastlines()
p2 = plt.contourf(lons, lats, oxy_o_dec, transform=ccrs.PlateCarree(), cmap=colmap2, levels=levs2, vmin=np.min(levs2),vmax=np.max(levs2), zorder=1, extend='both')
plt.contour(lons, lats, oxy_o_dec, transform=ccrs.PlateCarree(), colors=contcol, linewidths=contwid, levels=cont2, zorder=2)
'''

ax3 = plt.subplot(gs[12:20,0])
ax3.tick_params(labelsize=fstic)
plt.scatter(aou_m_dec, sol_m_dec, c=oxy_m_dec, cmap=colmap1, vmin=np.min(levs1), vmax=np.max(levs1), marker=mar[0], edgecolor='none', linewidths=0, alpha=0.5, zorder=zor[0])
c1 = plt.contour(aou_cont, sol_cont, oxy_cont, levels=np.arange(-100,101,5), colors=contcol, linewidths=contwid)
plt.plot((x1a,x2a),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1a,y2a), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1a,y2a); plt.xlim(x1a,x2a)
plt.ylabel('$\Delta$ O$_2^{sol}$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
plt.xlabel('$\Delta$ O$_2^{AOU}$ ($\mu$M decade$^{-1}$)', fontsize=fslab)


ax4 = plt.subplot(gs[12:20,1])
ax4.tick_params(labelsize=fstic, left=False, labelleft=False, right=True, labelright=True)
plt.scatter(aou_o_dec, sol_o_dec, c=oxy_o_dec, cmap=colmap2, vmin=np.min(levs2), vmax=np.max(levs2), marker=mar[0], edgecolor='none', linewidths=0, alpha=0.5, zorder=zor[0])
c2 = plt.contour(aou_cont, sol_cont, oxy_cont, levels=np.arange(-100,101,5), colors=contcol, linewidths=contwid)
plt.plot((x1b,x2b),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1b,y2b), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1b,y2b); plt.xlim(x1b,x2b)
plt.ylabel('$\Delta$ O$_2^{sol}$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
plt.xlabel('$\Delta$ O$_2^{AOU}$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
ax4.yaxis.set_label_position('right')


ax5 = plt.subplot(gs[24:32,0])
ax5.tick_params(labelsize=fstic)
plt.scatter(tou_m_dec, dis_m_dec, c=aou_m_dec, cmap=colmap1, vmin=np.min(levs1), vmax=np.max(levs1), marker=mar[0], edgecolor='none', linewidths=0, alpha=0.5, zorder=zor[0])
c3 = plt.contour(aou_cont, sol_cont, oxy_cont, levels=np.arange(-100,101,5), colors=contcol, linewidths=contwid)
plt.plot((x1b,x2b),(0,0), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.plot((0,0),(y1b,y2b), color=zerocol, linestyle=zerosty, linewidth=zerowid, zorder=zerozor)
plt.ylim(y1c,y2c); plt.xlim(x1c,x2c)
plt.ylabel('$\Delta$ O$_2^{dis}$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
plt.xlabel('$\Delta$ O$_2^{TOU}$ ($\mu$M decade$^{-1}$)', fontsize=fslab)

'''
ax6 = plt.subplot(gs[24:32,1])
ax6.spines['top'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax6.tick_params(labelsize=fstic, left=False, labelleft=False, right=True, labelright=True)
plt.plot(years_cmip6, obs*0, color=contcol, linewidth=contwid, linestyle=contsty, alpha=contalf, zorder=contzor)
plt.plot(years_cmip6, obs, color=fcol[0], linewidth=fwid[0], linestyle=fsty[0], alpha=falf[0], zorder=fzor[0], label=flab[0])
plt.plot(years_jra, jra, color=fcol[1], linewidth=fwid[1], linestyle=fsty[1], alpha=falf[1], zorder=fzor[1], label=flab[1])
plt.plot(years_cmip6[0:64], np.min(cmip6,axis=1), color=fcol[2], linewidth=fwid[2], linestyle=fsty[2], alpha=falf[2], zorder=fzor[2])
plt.plot(years_cmip6[0:64], np.max(cmip6,axis=1), color=fcol[2], linewidth=fwid[2], linestyle=fsty[2], alpha=falf[2], zorder=fzor[2])
plt.fill_between(years_cmip6[0:64], np.min(cmip6,axis=1), np.max(cmip6,axis=1), color=fcol[2], alpha=falf[2], label=flab[2])
plt.legend(frameon=False, loc='lower left', ncol=1, fontsize=11, labelspacing=0.25, columnspacing=0.5)
plt.xlim(1951,2015)
plt.ylim(-6,3)
plt.xlabel('Year', fontsize=fslab)
plt.ylabel('$\Delta$ O$_2$ ($\mu$M)', fontsize=fslab)
ax6.yaxis.set_label_position('right')
'''

fig.subplots_adjust(top=0.95, bottom=0.1, left=0.125, right=0.875, wspace=0.2, hspace=0.3)

xx = 0.5; yy = 1.05
plt.text(xx, yy, 'Hindcast simulation', fontsize=fstic, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'Hindcast simulation', fontsize=fstic, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'Observation-based product', fontsize=fstic, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'Hindcast simulation', fontsize=fstic, transform=ax5.transAxes, va='center', ha='center', fontweight='bold')


xx = 0.05; yy = 0.95
plt.text(xx, yy, 'a', fontsize=fslab+2, transform=ax1.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy, 'b', fontsize=fslab+2, transform=ax2.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy+0.1, 'c', fontsize=fslab+2, transform=ax3.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy+0.1, 'd', fontsize=fslab+2, transform=ax4.transAxes, va='center', ha='center', fontweight='bold')
plt.text(xx, yy+0.1, 'e', fontsize=fslab+2, transform=ax5.transAxes, va='center', ha='center', fontweight='bold')
#plt.text(xx, yy+0.1, 'f', fontsize=fslab+2, transform=ax6.transAxes, va='center', ha='center', fontweight='bold')


cbax1 = fig.add_axes([0.085, 0.68, 0.025, 0.22])
cbar1 = plt.colorbar(p1, cax=cbax1, orientation='vertical', ticks=levs1[::3])
cbax1.set_ylabel('$\Delta$ O$_2$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
cbax1.tick_params(labelsize=fstic, left=True, labelleft=True, right=False, labelright=False)
cbax1.yaxis.set_label_position('left')

'''
cbax2 = fig.add_axes([0.89, 0.68, 0.025, 0.22])
cbar2 = plt.colorbar(p2, cax=cbax2, orientation='vertical', ticks=levs2[::3])
cbax2.set_ylabel('$\Delta$ O$_2$ ($\mu$M decade$^{-1}$)', fontsize=fslab)
cbax2.tick_params(labelsize=fstic)
'''



#%% save figure

plt.clabel(c1, manual=True, fmt='%i', fontsize=fstic)
plt.clabel(c2, manual=True, fmt='%i', fontsize=fstic)
plt.clabel(c3, manual=True, fmt='%i', fontsize=fstic)
#plt.clabel(c4, manual=True, fmt='%.1f', fontsize=fstic)


os.chdir("C://Users//pearseb/Dropbox//PostDoc//my articles//historical model-data deoxygenation//final_figures")
fig.savefig('fig-fig1.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-fig1.eps', dpi=300)
fig.savefig('fig-fig1_trans.png', dpi=300, bbox_inches='tight', transparent=True)

