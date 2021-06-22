# -*- coding: utf-8 -*-
"""
Spyder Editor


make budgets figure

"""

#%% imports

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sb
sb.set(style='ticks')


#%% get data

os.chdir("C://Users/pearseb/Dropbox/PostDoc/my articles/historical model-data deoxygenation/data_for_figures")

data = nc.Dataset('figure_budget_regions.nc','r')

ben_tot = data.variables['BENGUELA_TOT_O2_SHFAVE_ZINT'][...]
ben_bio = data.variables['BENGUELA_BIO_O2_SHFAVE_ZINT'][...]
ben_phy = data.variables['BENGUELA_PHY_O2_SHFAVE_ZINT'][...]
ben_adv = data.variables['BENGUELA_ADV_O2_SHFAVE_ZINT'][...]
ben_hor = data.variables['BENGUELA_HOR_O2_SHFAVE_ZINT'][...]
ben_ver = data.variables['BENGUELA_VER_O2_SHFAVE_ZINT'][...]

can_tot = data.variables['CANARYATL_TOT_O2_SHFAVE_ZINT'][...]
can_bio = data.variables['CANARYATL_BIO_O2_SHFAVE_ZINT'][...]
can_phy = data.variables['CANARYATL_PHY_O2_SHFAVE_ZINT'][...]
can_adv = data.variables['CANARYATL_ADV_O2_SHFAVE_ZINT'][...]
can_hor = data.variables['CANARYATL_HOR_O2_SHFAVE_ZINT'][...]
can_ver = data.variables['CANARYATL_VER_O2_SHFAVE_ZINT'][...]

per_tot = data.variables['PERUCHILE_TOT_O2_SHFAVE_ZINT'][...]
per_bio = data.variables['PERUCHILE_BIO_O2_SHFAVE_ZINT'][...]
per_phy = data.variables['PERUCHILE_PHY_O2_SHFAVE_ZINT'][...]
per_adv = data.variables['PERUCHILE_ADV_O2_SHFAVE_ZINT'][...]
per_hor = data.variables['PERUCHILE_HOR_O2_SHFAVE_ZINT'][...]
per_ver = data.variables['PERUCHILE_VER_O2_SHFAVE_ZINT'][...]

cal_tot = data.variables['CALIFORNIA_TOT_O2_SHFAVE_ZINT'][...]
cal_bio = data.variables['CALIFORNIA_BIO_O2_SHFAVE_ZINT'][...]
cal_phy = data.variables['CALIFORNIA_PHY_O2_SHFAVE_ZINT'][...]
cal_adv = data.variables['CALIFORNIA_ADV_O2_SHFAVE_ZINT'][...]
cal_hor = data.variables['CALIFORNIA_HOR_O2_SHFAVE_ZINT'][...]
cal_ver = data.variables['CALIFORNIA_VER_O2_SHFAVE_ZINT'][...]

soc_tot = data.variables['SOUTHERNOCEAN_TOT_O2_SHFAVE_ZINT'][...]
soc_bio = data.variables['SOUTHERNOCEAN_BIO_O2_SHFAVE_ZINT'][...]
soc_phy = data.variables['SOUTHERNOCEAN_PHY_O2_SHFAVE_ZINT'][...]
soc_adv = data.variables['SOUTHERNOCEAN_ADV_O2_SHFAVE_ZINT'][...]
soc_hor = data.variables['SOUTHERNOCEAN_HOR_O2_SHFAVE_ZINT'][...]
soc_ver = data.variables['SOUTHERNOCEAN_VER_O2_SHFAVE_ZINT'][...]

ona_tot = data.variables['OLIGONA_TOT_O2_SHFAVE_ZINT'][...]
ona_bio = data.variables['OLIGONA_BIO_O2_SHFAVE_ZINT'][...]
ona_phy = data.variables['OLIGONA_PHY_O2_SHFAVE_ZINT'][...]
ona_adv = data.variables['OLIGONA_ADV_O2_SHFAVE_ZINT'][...]
ona_hor = data.variables['OLIGONA_HOR_O2_SHFAVE_ZINT'][...]
ona_ver = data.variables['OLIGONA_VER_O2_SHFAVE_ZINT'][...]

onp_tot = data.variables['OLIGONP_TOT_O2_SHFAVE_ZINT'][...]
onp_bio = data.variables['OLIGONP_BIO_O2_SHFAVE_ZINT'][...]
onp_phy = data.variables['OLIGONP_PHY_O2_SHFAVE_ZINT'][...]
onp_adv = data.variables['OLIGONP_ADV_O2_SHFAVE_ZINT'][...]
onp_hor = data.variables['OLIGONP_HOR_O2_SHFAVE_ZINT'][...]
onp_ver = data.variables['OLIGONP_VER_O2_SHFAVE_ZINT'][...]

ben = np.array([ben_tot, ben_bio, ben_phy, ben_adv, ben_hor, ben_ver])
can = np.array([can_tot, can_bio, can_phy, can_adv, can_hor, can_ver])
per = np.array([per_tot, per_bio, per_phy, per_adv, per_hor, per_ver])
cal = np.array([cal_tot, cal_bio, cal_phy, cal_adv, cal_hor, cal_ver])
soc = np.array([soc_tot, soc_bio, soc_phy, soc_adv, soc_hor, soc_ver])
ona = np.array([ona_tot, ona_bio, ona_phy, ona_adv, ona_hor, ona_ver])
onp = np.array([onp_tot, onp_bio, onp_phy, onp_adv, onp_hor, onp_ver])


data.close()


#%% what are the relative changes (demand : supply)

print("Benguela bio:sup = ",abs(ben_bio) / abs(ben_phy))
print("Canary bio:sup = ",abs(can_bio) / abs(can_phy))
print("Humboldt bio:sup = ",abs(per_bio) / abs(per_phy))
print("California bio:sup = ",abs(cal_bio) / abs(cal_phy))


#%% make figure


fstic = 13
fslab = 15

col = ['k', 'royalblue', 'firebrick', 'firebrick', 'firebrick', 'firebrick']
alf = [0.8, 0.8, 0.8, 0.5, 0.5, 0.5]
edc  = ['k', 'k', 'k', 'k', 'k', 'k']
lwi = [1,1,1,0.5,0.5,0.5]
zor = [1,1,1,1,1,1]
wid= [1,1,1,1,1,1]
lab = ['$\Delta$O$_2$', '$\Delta$O$_2^{bio}$', '$\Delta$O$_2^{phy}$', '$\Delta$O$_2^{adv}$', '$\Delta$O$_2^{hdf}$', '$\Delta$O$_2^{vdf}$']
hat = ['', '', '', '...', '---', '///']

fig = plt.figure(figsize=(10,6), facecolor='w')
gs = GridSpec(1,1)


ax1 = plt.subplot(gs[0])
ax1.tick_params(labelsize=fstic, bottom=False, labelbottom=False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

plt.plot((0,26.5),(0,0), 'k--', linewidth=0.5, zorder=1)

for i in np.arange(0,6):
    plt.bar(i, ben[i], facecolor=col[i], alpha=alf[i], edgecolor=edc[i], linewidth=lwi[i], zorder=zor[i], width=wid[i], hatch=hat[i])
    plt.bar(i+7, can[i], facecolor=col[i], alpha=alf[i], edgecolor=edc[i], linewidth=lwi[i], zorder=zor[i], width=wid[i], hatch=hat[i])
    plt.bar(i+14, per[i], facecolor=col[i], alpha=alf[i], edgecolor=edc[i], linewidth=lwi[i], zorder=zor[i], width=wid[i], hatch=hat[i])
    plt.bar(i+21, cal[i], facecolor=col[i], alpha=alf[i], edgecolor=edc[i], linewidth=lwi[i], zorder=zor[i], width=wid[i], hatch=hat[i], label=lab[i])
    #plt.bar(i+28, soc[i], facecolor=col[i], alpha=alf[i], edgecolor=edc[i], linewidth=lwi[i], zorder=zor[i], width=wid[i], hatch=hat[i], label=lab[i])

plt.ylim(-60,60)
plt.xlim(-1,27)

plt.legend(frameon=False, ncol=6, loc='lower center')

plt.text(2,65, 'Benguela\nupwelling', transform=ax1.transData, va='center', ha='center', rotation=45, fontsize=fslab)
plt.text(9,65, 'Canary\nupwelling', transform=ax1.transData, va='center', ha='center', rotation=45, fontsize=fslab)
plt.text(16,65, 'Peru/Chile\nupwelling', transform=ax1.transData, va='center', ha='center', rotation=45, fontsize=fslab)
plt.text(23,65, 'California\nupwelling', transform=ax1.transData, va='center', ha='center', rotation=45, fontsize=fslab)
#plt.text(30,65, 'Southern\nOcean', transform=ax1.transData, va='center', ha='center', rotation=45, fontsize=fslab)

plt.ylabel('$\Delta$O$_2$ (mmol m$^{-2}$ year$^{-1}$)', fontsize=fslab)


#%% save figure

os.chdir("C://Users/pearseb/Dropbox/PostDoc/my articles/historical model-data deoxygenation/final_figures")
fig.savefig('fig-suppfig4.png', dpi=300, bbox_inches='tight')
fig.savefig('fig-suppfig4_trans.png', dpi=300, bbox_inches='tight', transparent=True)


