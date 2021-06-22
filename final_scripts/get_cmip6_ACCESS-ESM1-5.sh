#!/bin/bash

export WRK=/mnt/lustre/users/pearseb/NEMO_OUT/analysis_oxygen/cmip6
mkdir -p $WRK
cd $WRK

source /users/pearseb/purge.env

	# ACCESS-ESM1-5 - uses MOM5 ocean model
### oxygen
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Oyr/o2/gn/v20191115/o2_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1950-2014.nc
ncks -O -F -d time,2,65,1 o2_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1950-2014.nc o2_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
rm o2_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1950-2014.nc

### NPP
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/intpp/gn/v20191115/intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc intpp_${year}.nc
 cdo -O -f nc -yearavg intpp_${year}.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncrcat ${files} intpp_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} intpp_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc

### Temperature
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_195001-195912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_196001-196912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_197001-197912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_198001-198912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_199001-199912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_200001-200912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/thetao/gn/v20191115/thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls thetao_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_*.nc`
source /users/pearseb/purge.env
ncrcat -O ${files} thetao.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} thetao.nc thetao_${year}.nc
 ncra -O thetao_${year}.nc thetao_${year}.nc
done
rm ${files}
files=`ls thetao_????.nc`
source /users/pearseb/purge.env
ncrcat ${files} thetao_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} thetao.nc

wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_195001-195912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_196001-196912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_197001-197912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_198001-198912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_199001-199912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_200001-200912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/so/gn/v20191115/so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls so_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_*.nc`
source /users/pearseb/purge.env
ncrcat -O ${files} so.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} so.nc so_${year}.nc
 ncra -O so_${year}.nc so_${year}.nc
done
rm ${files}
files=`ls so_????.nc`
source /users/pearseb/purge.env
ncrcat ${files} so_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} so.nc 

wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_195001-195912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_196001-196912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_197001-197912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_198001-198912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_199001-199912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_200001-200912.nc
wget http://esgf.nci.org.au/thredds/fileServer/master/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/agessc/gn/v20191115/agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls agessc_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_*.nc`
source /users/pearseb/purge.env
ncrcat -O ${files} agessc.nc
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} agessc.nc agessc_${year}.nc
 ncra -O agessc_${year}.nc agessc_${year}.nc
done
rm ${files}
files=`ls agessc_????.nc`
source /users/pearseb/purge.env
ncrcat ${files} agessc_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} agessc.nc



### do some regridding
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/no_land_regrid.sh o2_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh so_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh thetao_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh intpp_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh agessc_Oyr_ACCESS-ESM1-5_historical_r1i1p1f1_gn_1951-2014.nc




