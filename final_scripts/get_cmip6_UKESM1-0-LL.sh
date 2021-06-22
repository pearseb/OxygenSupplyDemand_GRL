#!/bin/bash

export WRK=/mnt/lustre/users/pearseb/NEMO_OUT/analysis_oxygen/cmip6
mkdir -p $WRK
cd $WRK

source /users/pearseb/purge.env

	# UKESM1-0-LL - uses MRICOM4 ocean model
### oxygen
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/o2/gn/v20190627/o2_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/o2/gn/v20190627/o2_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc
files=`ls o2_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_*.nc`
ncrcat ${files} o2_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} o2_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc o2_${year}.nc
 ncra -O o2_${year}.nc o2_${year}.nc
done
source /users/pearseb/purge.env
files=`ls o2_????.nc`
ncrcat ${files} o2_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} o2_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

### NPP
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/intpp/gn/v20190627/intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc intpp_${year}.nc
 ncra -O intpp_${year}.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncrcat ${files} intpp_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} intpp_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

### Age
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/agessc/gn/v20190627/agessc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/agessc/gn/v20190627/agessc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc
files=`ls agessc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_*.nc`
ncrcat ${files} agessc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} agessc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc agessc_${year}.nc
 ncra -O agessc_${year}.nc agessc_${year}.nc
done
source /users/pearseb/purge.env
files=`ls agessc_????.nc`
ncrcat ${files} agessc_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} agessc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

### Temperature
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/thetao/gn/v20190627/thetao_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/thetao/gn/v20190627/thetao_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc
files=`ls thetao_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_*.nc`
ncrcat ${files} thetao_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} thetao_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc thetao_${year}.nc
 ncra -O thetao_${year}.nc thetao_${year}.nc
done
source /users/pearseb/purge.env
files=`ls thetao_????.nc`
ncrcat ${files} thetao_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} thetao_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

### Salinity
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/so/gn/v20190627/so_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/so/gn/v20190627/so_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc
files=`ls so_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_*.nc`
ncrcat ${files} so_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} so_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc so_${year}.nc
 ncra -O so_${year}.nc so_${year}.nc
done
source /users/pearseb/purge.env
files=`ls so_????.nc`
ncrcat ${files} so_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} so_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

### carbon export
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/expc/gn/v20190627/expc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/expc/gn/v20190627/expc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc
files=`ls expc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_*.nc`
ncrcat ${files} expc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} expc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc expc_${year}.nc
 ncra -O expc_${year}.nc expc_${year}.nc
done
source /users/pearseb/purge.env
files=`ls expc_????.nc`
ncrcat ${files} expc_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} expc_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

### vertical velocity
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/wo/gn/v20190627/wo_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc
wget http://esgf-data3.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/wo/gn/v20190627/wo_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc
files=`ls wo_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_*.nc`
ncrcat ${files} wo_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} wo_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc wo_${year}.nc
 ncra -O wo_${year}.nc wo_${year}.nc
done
source /users/pearseb/purge.env
files=`ls wo_????.nc`
ncrcat ${files} wo_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} wo_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc



### do some regridding
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/no_land_regrid.sh o2_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh so_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh thetao_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh intpp_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh agessc_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh expc_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh wo_Oyr_UKESM1-0-LL_historical_r1i1p1f2_gn_1951-2014.nc




