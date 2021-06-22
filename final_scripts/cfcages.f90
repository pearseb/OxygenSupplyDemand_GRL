program CFCage
  use netcdf

  implicit none
  
 
  character(len=40), parameter :: file_cfc11 = "ETOPO_JRA55_pic_1m_pCFC11_T_1958-2019.nc"
  character(len=40), parameter :: file_cfc12 = "ETOPO_JRA55_pic_1m_pCFC12_T_1958-2019.nc"
  character(len=38), parameter :: file_sf6 = "ETOPO_JRA55_pic_1m_pSF6_T_1958-2019.nc"
  integer :: ncid, varid

  ! dimensions
  integer, parameter :: ndims = 5
  integer, parameter :: nyears = 62
  integer, parameter :: nmonths = 12
  integer, parameter :: ndeps = 31
  integer, parameter :: nlats = 180
  integer, parameter :: nlons = 360
  ! loop indices
  integer :: i,j,k,l,m
  real :: years(nyears), months(nmonths), deps(ndeps), lats(nlats), lons(nlons)
  ! variables to be read in
  real :: pcfc11(nlons,nlats,ndeps,nmonths,nyears)
  real :: pcfc12(nlons,nlats,ndeps,nmonths,nyears)
  real :: psf6(nlons,nlats,ndeps,nmonths,nyears)
  real :: atmdata(89,7)
  ! variables to be calculated
  real :: cfc11_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: cfc12_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: sf6_age(nlons,nlats,ndeps,nmonths,nyears)
  ! temporary variables required during calculation
  real :: z_val
  real :: i_val(1), i2_val(1)
  real :: n_val(1)
  real :: y1_val(1), y2_val(1), v1_val(1), v2_val(1), y_val(1)
  ! for saving netcdf files with age data
  integer :: yeardim_id, monthdim_id, depdim_id, latdim_id, londim_id
  integer :: yearvar_id, monthvar_id, depvar_id, latvar_id, lonvar_id
  integer :: cfc11var_id, cfc12var_id, sf6var_id
  integer :: dimids(5)
  

  !!! GET NETCDF DATA !!!

  ! open file 
  call check( nf90_open(file_cfc11, nf90_nowrite, ncid) )
  ! get variable id
  call check( nf90_inq_varid(ncid, "PCFC11", varid) )
  ! read data
  call check( nf90_get_var(ncid, varid, pcfc11(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  !close file
  call check (nf90_close(ncid) )
  
  ! open file 
  call check( nf90_open(file_cfc12, nf90_nowrite, ncid) )
  ! get variable id
  call check( nf90_inq_varid(ncid, "PCFC12", varid) )
  ! read data
  call check( nf90_get_var(ncid, varid, pcfc12(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  !close file
  call check (nf90_close(ncid) )
  
  ! open file 
  call check( nf90_open(file_sf6, nf90_nowrite, ncid) )
  ! get variable id
  call check( nf90_inq_varid(ncid, "PSF6", varid) )
  ! read data
  call check( nf90_get_var(ncid, varid, psf6(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  !close file
  call check (nf90_close(ncid) )


  ! define years
  do m = 1,nyears
   years(m) = 1958.5 + (m-1)
  enddo
  ! define months
  do l = 1,nmonths
   months(l) = 1 + (l-1)
  enddo
  ! define depths
  deps = (/ 5,15,25,35,45,55,65,75,85,95,106,117,129,142,159,182,217,272,364,  &
            512,732,1033,1406,1831,2290,2768,3257,3752,4250,4750,5250 /)
  ! define lats
  do j = 1,nlats
   lats(j) = -89.5 + (j-1)
  enddo
  ! define lons
  do i = 1,nlons
   lons(i) = -89.5 + (i-1)
  enddo
 
  !!! read in atmospheric CFC data !!!
  open(unit=10, file='CFCs_CDIAC.dat', form='formatted', status='old', action='read')
  do i = 1,89
    read(10, fmt=*) atmdata(i,:)
  enddo 
  close(10)

  print*, " *** "
  print*, " Doing CFC-11 "
  !!! Conduct calculation of CFC-11 ages !!!
  do i = 1,nlons
  print*, " i ", i
   do j = 1,nlats
    do k = 1,ndeps
     do l = 1,nmonths
      do m = 1,nyears

        ! make sure all negative values are zero
        if (pcfc11(i,j,k,l,m).lt.0.0) pcfc11(i,j,k,l,m) = 0.0     
      
        if (years(m).lt.1995) then  ! if year is less than 1995 do cfc11
           
          z_val = pcfc11(i,j,k,l,m)

          if (lats(j).gt.0) then   ! Northern Hemisphere data
            if (z_val.ge.maxval(atmdata(:,2))) then
              cfc11_age(i,j,k,l,m) = 0.0
            else
              i_val = minloc(abs(atmdata(:,2) - z_val))
              n_val = atmdata(i_val,2) - z_val
              i2_val = i_val - sign(1.0,n_val(1))
              
              y1_val = atmdata(i_val,1)
              y2_val = atmdata(i2_val,1)
              v1_val = atmdata(i_val,2)
              v2_val = atmdata(i2_val,2)
            endif
          else  ! Southern Hemisphere data
            if (z_val.ge.maxval(atmdata(:,3))) then
              cfc11_age(i,j,k,l,m) = 0.0
            else
              i_val = minloc(abs(atmdata(:,3) - z_val))
              n_val = atmdata(i_val,3) - z_val
              i2_val = i_val - sign(1.0,n_val(1))
              
              y1_val = atmdata(i_val,1)
              y2_val = atmdata(i2_val,1)
              v1_val = atmdata(i_val,3)
              v2_val = atmdata(i2_val,3)
            endif
          endif

          y_val = y1_val + (n_val / (v1_val - v2_val)) &
                    * (y2_val - y1_val)
          cfc11_age(i,j,k,l,m) = years(m) - y_val(1)
        
        else  ! years >= 1995

          cfc11_age(i,j,k,l,m) = -9999.0

        endif

      enddo  ! m
     enddo  ! l
    enddo  ! k
   enddo  ! j
  enddo  ! i


  print*, " *** "
  print*, " Doing CFC-12 "
  !!! Conduct calculation of CFC-12 ages !!!
  do i = 1,nlons
  print*, " i ", i
   do j = 1,nlats
    do k = 1,ndeps
     do l = 1,nmonths
      do m = 1,nyears

        ! make sure all negative values are zero
        if (pcfc12(i,j,k,l,m).lt.0.0) pcfc12(i,j,k,l,m) = 0.0     
      
        if (years(m).le.2004) then  ! if year is less than 2005 do cfc12
           
          z_val = pcfc12(i,j,k,l,m)

          if (lats(j).gt.0) then   ! Northern Hemisphere data
            if (z_val.ge.maxval(atmdata(:,4))) then
              cfc12_age(i,j,k,l,m) = 0.0
            else
              i_val = minloc(abs(atmdata(:,4) - z_val))
              n_val = atmdata(i_val,4) - z_val
              i2_val = i_val - sign(1.0,n_val(1))
              
              y1_val = atmdata(i_val,1)
              y2_val = atmdata(i2_val,1)
              v1_val = atmdata(i_val,4)
              v2_val = atmdata(i2_val,4)
            endif
          else  ! Southern Hemisphere data
            if (z_val.ge.maxval(atmdata(:,5))) then
              cfc12_age(i,j,k,l,m) = 0.0
            else
              i_val = minloc(abs(atmdata(:,5) - z_val))
              n_val = atmdata(i_val,5) - z_val
              i2_val = i_val - sign(1.0,n_val(1))
              
              y1_val = atmdata(i_val,1)
              y2_val = atmdata(i2_val,1)
              v1_val = atmdata(i_val,5)
              v2_val = atmdata(i2_val,5)
            endif
          endif

          y_val = y1_val + (n_val / (v1_val - v2_val)) &
                    * (y2_val - y1_val)
          cfc12_age(i,j,k,l,m) = years(m) - y_val(1)
        
        else  ! years > 2004

          cfc12_age(i,j,k,l,m) = -9999.0

        endif

      enddo  ! m
     enddo  ! l
    enddo  ! k
   enddo  ! j
  enddo  ! i


  print*, " *** "
  print*, " Doing SF6 "
  !!! Conduct calculation of SF6 ages !!!
  do i = 1,nlons
  print*, " i ", i
   do j = 1,nlats
    do k = 1,ndeps
     do l = 1,nmonths
      do m = 1,nyears

        ! make sure all negative values are zero
        if (psf6(i,j,k,l,m).lt.0.0) psf6(i,j,k,l,m) = 0.0     
      
        if (years(m).ge.1990) then  ! if year gt 1989 do sf6
           
          z_val = psf6(i,j,k,l,m)

          if (lats(j).gt.0) then   ! Northern Hemisphere data
            if (z_val.ge.maxval(atmdata(:,6))) then
              sf6_age(i,j,k,l,m) = 0.0
            else
              i_val = minloc(abs(atmdata(:,6) - z_val))
              n_val = atmdata(i_val,6) - z_val
              i2_val = i_val - sign(1.0,n_val(1))
              
              y1_val = atmdata(i_val,1)
              y2_val = atmdata(i2_val,1)
              v1_val = atmdata(i_val,6)
              v2_val = atmdata(i2_val,6)
            endif
          else  ! Southern Hemisphere data
            if (z_val.ge.maxval(atmdata(:,7))) then
              sf6_age(i,j,k,l,m) = 0.0
            else
              i_val = minloc(abs(atmdata(:,7) - z_val))
              n_val = atmdata(i_val,7) - z_val
              i2_val = i_val - sign(1.0,n_val(1))
              
              y1_val = atmdata(i_val,1)
              y2_val = atmdata(i2_val,1)
              v1_val = atmdata(i_val,7)
              v2_val = atmdata(i2_val,7)
            endif
          endif

          y_val = y1_val + (n_val / (v1_val - v2_val)) &
                    * (y2_val - y1_val)
          sf6_age(i,j,k,l,m) = years(m) - y_val(1)
        
        else  ! years < 1990

          sf6_age(i,j,k,l,m) = -9999.0

        endif

      enddo  ! m
     enddo  ! l
    enddo  ! k
   enddo  ! j
  enddo  ! i

  print*, " *** "
  print*, " Saving data to new netcdf files "
  !!! Save the data
  call check( nf90_create("ETOPO_JRA55_pic_1m_CFC11age_T_1958-2019.nc", nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'RECORD', nyears, yeardim_id) )
  call check( nf90_def_dim(ncid, 'TIME_COUNTER', nmonths, monthdim_id) )
  call check( nf90_def_dim(ncid, 'DEPTHT', ndeps, depdim_id) )
  call check( nf90_def_dim(ncid, 'ETOPO60Y', nlats, latdim_id) )
  call check( nf90_def_dim(ncid, 'ETOPO60X', nlons, londim_id) )
  call check( nf90_def_var(ncid, 'year', nf90_real, yeardim_id, yearvar_id) )
  call check( nf90_def_var(ncid, 'month', nf90_real, monthdim_id, monthvar_id) )
  call check( nf90_def_var(ncid, 'depth', nf90_real, depdim_id, depvar_id) )
  call check( nf90_def_var(ncid, 'latitude', nf90_real, latdim_id, latvar_id) )
  call check( nf90_def_var(ncid, 'longitude', nf90_real, londim_id, lonvar_id) )
  dimids = (/ londim_id, latdim_id, depdim_id, monthdim_id, yeardim_id /)
  call check( nf90_def_var(ncid, 'CFC11age', nf90_real, dimids, cfc11var_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, cfc11var_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, cfc11var_id, cfc11_age, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  call check( nf90_create("ETOPO_JRA55_pic_1m_CFC12age_T_1958-2019.nc", nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'RECORD', nyears, yeardim_id) )
  call check( nf90_def_dim(ncid, 'TIME_COUNTER', nmonths, monthdim_id) )
  call check( nf90_def_dim(ncid, 'DEPTHT', ndeps, depdim_id) )
  call check( nf90_def_dim(ncid, 'ETOPO60Y', nlats, latdim_id) )
  call check( nf90_def_dim(ncid, 'ETOPO60X', nlons, londim_id) )
  call check( nf90_def_var(ncid, 'year', nf90_real, yeardim_id, yearvar_id) )
  call check( nf90_def_var(ncid, 'month', nf90_real, monthdim_id, monthvar_id) )
  call check( nf90_def_var(ncid, 'depth', nf90_real, depdim_id, depvar_id) )
  call check( nf90_def_var(ncid, 'latitude', nf90_real, latdim_id, latvar_id) )
  call check( nf90_def_var(ncid, 'longitude', nf90_real, londim_id, lonvar_id) )
  dimids = (/ londim_id, latdim_id, depdim_id, monthdim_id, yeardim_id /)
  call check( nf90_def_var(ncid, 'CFC12age', nf90_real, dimids, cfc12var_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, cfc12var_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, cfc12var_id, cfc12_age, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  call check( nf90_create("ETOPO_JRA55_pic_1m_SF6age_T_1958-2019.nc", nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'RECORD', nyears, yeardim_id) )
  call check( nf90_def_dim(ncid, 'TIME_COUNTER', nmonths, monthdim_id) )
  call check( nf90_def_dim(ncid, 'DEPTHT', ndeps, depdim_id) )
  call check( nf90_def_dim(ncid, 'ETOPO60Y', nlats, latdim_id) )
  call check( nf90_def_dim(ncid, 'ETOPO60X', nlons, londim_id) )
  call check( nf90_def_var(ncid, 'year', nf90_real, yeardim_id, yearvar_id) )
  call check( nf90_def_var(ncid, 'month', nf90_real, monthdim_id, monthvar_id) )
  call check( nf90_def_var(ncid, 'depth', nf90_real, depdim_id, depvar_id) )
  call check( nf90_def_var(ncid, 'latitude', nf90_real, latdim_id, latvar_id) )
  call check( nf90_def_var(ncid, 'longitude', nf90_real, londim_id, lonvar_id) )
  dimids = (/ londim_id, latdim_id, depdim_id, monthdim_id, yeardim_id /)
  call check( nf90_def_var(ncid, 'SF6age', nf90_real, dimids, sf6var_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, sf6var_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, sf6var_id, sf6_age, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )


contains

  subroutine check(status)
    integer, intent(in) :: status
    
    if (status /= nf90_noerr) then
     print*, trim(nf90_strerror(status))
     stop "Stopped"
    endif
  end subroutine check

end program CFCage
