#!/bin/bash

### download and process model output (includes regridding)
bash get_cmip6_ACCESS-ESM1-5.sh
bash get_cmip6_CanESM5.sh
#bash get_cmip6_CESM2.sh
#bash get_cmip6_CESM2-WACCM.sh
#bash get_cmip6_CNRM-ESM2-1.sh
bash get_cmip6_IPSL-CM6A-LR.sh
bash get_cmip6_MIROC-ES2L.sh
bash get_cmip6_MPI-ESM1-2-LR.sh
bash get_cmip6_MRI-ESM2-0.sh
bash get_cmip6_UKESM1-0-LL.sh


### calculate o2sat from thetao and so
ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl ETOPO_thetao_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc ETOPO_so_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc ETOPO_o2sat_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
#ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl  
ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl ETOPO_thetao_Oyr_IPSL-CM6A-LR_historical_r1i1p1f1_gn_1951-2014.nc ETOPO_so_Oyr_IPSL-CM6A-LR_historical_r1i1p1f1_gn_1951-2014.nc ETOPO_o2sat_Oyr_IPSL-CM6A-LR_historical_r1i1p1f1_gn_1951-2014.nc
ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl ETOPO_thetao_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc ETOPO_so_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc ETOPO_o2sat_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc 
ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl ETOPO_thetao_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc ETOPO_so_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc ETOPO_o2sat_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl ETOPO_thetao_Oyr_MRI-ESM2-0_historical_r1i2p1f1_gn_1951-2014.nc ETOPO_so_Oyr_MRI-ESM2-0_historical_r1i2p1f1_gn_1951-2014.nc ETOPO_o2sat_Oyr_MRI-ESM2-0_historical_r1i2p1f1_gn_1951-2014.nc
ferret -script /users/pearseb/analysis_oxygen/cmip6_o2sat.jnl ETOPO_thetao_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc ETOPO_so_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc ETOPO_o2sat_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc


### calculate linear global trends
ferret -script /users/pearseb/analysis_oxygen/cmip6_multimodel_trends.jnl

### calculate historical trends (2005-2014 minus 1975-1984), correlations with AOU and slopes and save to file
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_ACCESS-ESM1-5.jnl
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_CanESM5.jnl
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_IPSL-CM6A-LR.jnl
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_MIROC-ES2L.jnl
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_MPI-ESM1-2-LR.jnl
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_MRI-ESM2-0.jnl
ferret -script /users/pearseb/analysis_oxygen/cmip6_trends_correlations_UKESM1-0-LL.jnl

### view differences in AOU, vertical velocity, NPP an Cexp
ferret -script /users/pearseb/analysis_oxygen/cmip6_ebus.jnl

### save vertically integrated NPP and carbon export through 200 metres
