
rm cfcages.exe
rm cfcages_obs.exe
rm combine_data.exe

ifort -mcmodel=medium -I${NETCDFDIR}/include -I${NETCDFFORTRANDIR}/include -L${NETCDFFORTRANDIR}/lib -L${NETCDFDIR}/lib -lnetcdf -lnetcdff cfcages.f90 -o cfcages.exe
ifort -mcmodel=medium -I${NETCDFDIR}/include -I${NETCDFFORTRANDIR}/include -L${NETCDFFORTRANDIR}/lib -L${NETCDFDIR}/lib -lnetcdf -lnetcdff cfcages_obs.f90 -o cfcages_obs.exe
ifort -mcmodel=medium -I${NETCDFDIR}/include -I${NETCDFFORTRANDIR}/include -L${NETCDFFORTRANDIR}/lib -L${NETCDFDIR}/lib -lnetcdf -lnetcdff combine_data.f90 -o combine_data.exe

chmod 755 cfcages.exe
chmod 755 cfcages_obs.exe
chmod 755 combine_data.exe
