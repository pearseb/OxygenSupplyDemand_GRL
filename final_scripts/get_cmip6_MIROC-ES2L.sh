#!/bin/bash

export WRK=/mnt/lustre/users/pearseb/NEMO_OUT/analysis_oxygen/cmip6
mkdir -p $WRK
cd $WRK

source /users/pearseb/purge.env

	# MIROC-ES2L - uses COCO ocean model
### oxygen
wget http://esgf-data2.diasjp.net/thredds/fileServer/esg_dataroot/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Oyr/o2/gn/v20190823/o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1850-2014.nc
ncks -O -F -d time,102,165,1 o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1850-2014.nc o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
rm o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1850-2014.nc

### NPP
wget http://esgf-data2.diasjp.net/thredds/fileServer/esg_dataroot/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/intpp/gn/v20191129/intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc intpp_${year}.nc
 ncra -O intpp_${year}.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncrcat ${files} intpp_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc


### Age
wget http://esgf-data2.diasjp.net/thredds/fileServer/esg_dataroot/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/agessc/gn/v20190823/agessc_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc 
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} agessc_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc agessc_${year}.nc
 ncra -O agessc_${year}.nc agessc_${year}.nc
done
source /users/pearseb/purge.env
files=`ls agessc_????.nc`
ncrcat ${files} agessc_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} agessc_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc

### Temperature
wget http://esgf-data2.diasjp.net/thredds/fileServer/esg_dataroot/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/thetao/gn/v20190823/thetao_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc 
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} thetao_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc thetao_${year}.nc
 ncra -O thetao_${year}.nc thetao_${year}.nc
done
source /users/pearseb/purge.env
files=`ls thetao_????.nc`
ncrcat ${files} thetao_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} thetao_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc

### Salinity
wget http://esgf-data2.diasjp.net/thredds/fileServer/esg_dataroot/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/so/gn/v20190823/so_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc 
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} so_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc so_${year}.nc
 ncra -O so_${year}.nc so_${year}.nc
done
source /users/pearseb/purge.env
files=`ls so_????.nc`
ncrcat ${files} so_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
rm ${files} so_Omon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc


ncks -O -C -x -v depth,depth_c,eta,nsigma,sigma,sigma_bnds,zlev,zlev_bnds o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
ncks -O -C -x -v depth,depth_c,eta,nsigma,sigma,sigma_bnds,zlev,zlev_bnds agessc_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc agessc_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
ncks -O -C -x -v depth,depth_c,eta,nsigma,sigma,sigma_bnds,zlev,zlev_bnds thetao_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc thetao_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
ncks -O -C -x -v depth,depth_c,eta,nsigma,sigma,sigma_bnds,zlev,zlev_bnds so_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc so_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc

### do some regridding
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/no_land_regrid.sh o2_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh so_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh thetao_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh intpp_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh agessc_Oyr_MIROC-ES2L_historical_r1i1p1f2_gn_1951-2014.nc




