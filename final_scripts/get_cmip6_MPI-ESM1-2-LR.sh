#!/bin/bash

export WRK=/mnt/lustre/users/pearseb/NEMO_OUT/analysis_oxygen/cmip6
mkdir -p $WRK
cd $WRK

source /users/pearseb/purge.env

	# MPI-ESM1-2-LR - uses MPI ocean model
### oxygen
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/o2/gn/v20190710/o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-1969.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/o2/gn/v20190710/o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1970-1989.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/o2/gn/v20190710/o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1990-2009.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/o2/gn/v20190710/o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_2010-2014.nc
files=`ls o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat -O ${files} o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-2014.nc
ncks -O -F -d time,2,65,1 o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-2014.nc o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-2014.nc


### NPP
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-196912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_197001-198912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/intpp/gn/v20190710/intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat ${files} intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc intpp_${year}.nc
 ncra -O intpp_${year}.nc intpp_${year}.nc
done
source /users/pearseb/purge.env
files=`ls intpp_????.nc`
ncrcat ${files} intpp_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} intpp_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc


### Age
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/agessc/gn/v20190710/agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-196912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/agessc/gn/v20190710/agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_197001-198912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/agessc/gn/v20190710/agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/agessc/gn/v20190710/agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat ${files} agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc agessc_${year}.nc
 ncra -O agessc_${year}.nc agessc_${year}.nc
done
source /users/pearseb/purge.env
files=`ls agessc_????.nc`
ncrcat ${files} agessc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} agessc_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc

### Temperature
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/thetao/gn/v20190710/thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-196912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/thetao/gn/v20190710/thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_197001-198912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/thetao/gn/v20190710/thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/thetao/gn/v20190710/thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat ${files} thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc thetao_${year}.nc
 ncra -O thetao_${year}.nc thetao_${year}.nc
done
source /users/pearseb/purge.env
files=`ls thetao_????.nc`
ncrcat ${files} thetao_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} thetao_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc

### Salinity
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/so/gn/v20190710/so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-196912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/so/gn/v20190710/so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_197001-198912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/so/gn/v20190710/so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/so/gn/v20190710/so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat ${files} so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc so_${year}.nc
 ncra -O so_${year}.nc so_${year}.nc
done
source /users/pearseb/purge.env
files=`ls so_????.nc`
ncrcat ${files} so_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} so_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc

### carbon export
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/expc/gn/v20190710/expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-1969.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/expc/gn/v20190710/expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1970-1989.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/expc/gn/v20190710/expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1990-2009.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Oyr/expc/gn/v20190710/expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_2010-2014.nc
files=`ls expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat -O ${files} expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-2014.nc
ncks -O -F -d time,2,65,1 expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-2014.nc expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1950-2014.nc

### vertical velocity
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/wo/gn/v20190710/wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-196912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/wo/gn/v20190710/wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_197001-198912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/wo/gn/v20190710/wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_199001-200912.nc
wget http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/wo/gn/v20190710/wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_201001-201412.nc
files=`ls wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_*.nc`
ncrcat ${files} wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc
rm ${files}
source /users/pearseb/load_cdo.env
for year in $(seq 1951 1 2014); do
 cdo -O -f nc -selyear,${year} wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc wo_${year}.nc
 ncra -O wo_${year}.nc wo_${year}.nc
done
source /users/pearseb/purge.env
files=`ls wo_????.nc`
ncrcat ${files} wo_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
rm ${files} wo_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_195001-201412.nc



### do some regridding
source /users/pearseb/load_cdo.env
/users/pearseb/regridfiles/no_land_regrid.sh o2_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh so_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh thetao_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh intpp_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh agessc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh expc_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc
/users/pearseb/regridfiles/no_land_regrid.sh wo_Oyr_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_1951-2014.nc




