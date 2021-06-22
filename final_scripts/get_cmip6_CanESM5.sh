#!/bin/bash

export WRK=/mnt/lustre/users/pearseb/NEMO_OUT/analysis_oxygen/cmip6
mkdir -p $WRK
cd $WRK

source /users/pearseb/purge.env

	# CanESM5 - NEMO ocean model v3.4.1
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Oyr/o2/gn/v20190429/o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -O -F -d time,102,165,1 o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc
rm o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc

wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgA_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Oyr/o2sat/gn/v20190429/o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -O -F -d time,102,165,1 o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc
rm o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc

wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgA_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Oyr/pp/gn/v20190429/pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -O -F -d time,102,165,1 pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc
rm pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc

wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_195101-196012.nc
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_196101-197012.nc
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_197101-198012.nc
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_198101-199012.nc
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_199101-200012.nc
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_200101-201012.nc
wget http://crd-esgf-drc.ec.gc.ca/thredds/fileServer/esgC_dataroot/AR6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p2f1/Omon/agessc/gn/v20190429/agessc_Omon_CanESM5_historical_r1i1p2f1_gn_201101-201412.nc
files=`ls agessc_Omon_CanESM5_historical_r1i1p2f1_gn_*.nc`
source /users/pearseb/purge.env
ncrcat -O ${files} tmp.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} tmp.nc tmp_${year}.nc
done
source /users/pearseb/purge.env
for year in $(seq 1951 1 2014); do
 ncra -O tmp_${year}.nc tmp_${year}.nc
done
rm ${files}
files=`ls tmp_????.nc`
ncrcat ${files} agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc
rm ${files} tmp.nc agessc_Omon_CanESM5_historical_r1i1p2f1_gn_??????-??????.nc


ncks -A -v longitude agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -A -v latitude agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -A -v longitude agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -A -v latitude agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -A -v longitude agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc
ncks -A -v latitude agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc

### do some regridding
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/no_land_regrid.sh o2_Oyr_CanESM5_historical_r1i1p2f1_gn_1850-2014.nc 
/users/pearseb/regridfiles/no_land_regrid.sh o2sat_Oyr_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh pp_Oyr_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh agessc_Omon_CanESM5_historical_r1i1p2f1_gn_1951-2014.nc

