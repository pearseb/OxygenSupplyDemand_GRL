program CFCage
  use netcdf

  implicit none
  
 
  character(len=37), parameter :: fmod_wma = "ETOPO_JRA55_pic_1m_Age_T_1972-2019.nc"
  character(len=40), parameter :: fmod_oxy = "ETOPO_JRA55_pic_1m_oxygen_T_1972-2019.nc"
  character(len=45), parameter :: fmod_tem = "ETOPO_JRA55_pic_1m_temperature_T_1972-2019.nc"
  character(len=42), parameter :: fmod_sal = "ETOPO_JRA55_pic_1m_salinity_T_1972-2019.nc"
  character(len=42), parameter :: fmod_cfc11 = "ETOPO_JRA55_pic_1m_CFC11age_T_1972-2019.nc"
  character(len=42), parameter :: fmod_cfc12 = "ETOPO_JRA55_pic_1m_CFC12age_T_1972-2019.nc"
  character(len=40), parameter :: fmod_sf6 = "ETOPO_JRA55_pic_1m_SF6age_T_1972-2019.nc"
  character(len=43), parameter :: fmod_rem = "ETOPO_JRA55_pic_1m_reminnitr_T_1972-2019.nc"
  character(len=34), parameter :: fobs_oxy = "GLODAPv2.2020_oxygen_1x1degrees.nc"
  character(len=39), parameter :: fobs_tem = "GLODAPv2.2020_temperature_1x1degrees.nc"
  character(len=36), parameter :: fobs_sal = "GLODAPv2.2020_salinity_1x1degrees.nc"
  character(len=25), parameter :: fobs_cfc11 = "GLODAPv2.2020_CFC11age.nc"
  character(len=25), parameter :: fobs_cfc12 = "GLODAPv2.2020_CFC12age.nc"
  character(len=23), parameter :: fobs_sf6 = "GLODAPv2.2020_SF6age.nc"
  integer :: ncid, varid

  ! dimensions
  integer, parameter :: ndims = 5
  integer, parameter :: nyears = 48
  integer, parameter :: nmonths = 12
  integer, parameter :: ndeps = 31
  integer, parameter :: nlats = 180
  integer, parameter :: nlons = 360
  ! loop indices
  integer :: i,j,k,l,m, cnt,cnt2
  real :: years(nyears), months(nmonths), deps(ndeps), lats(nlats), lons(nlons)
  ! variables to be read in
  real :: mod_wma(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_cfc11_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_cfc12_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_sf6_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_oxy(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_ao2(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_tem(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_sal(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_rem(nlons,nlats,ndeps,nmonths,nyears)
  real :: mod_nit(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_cfc11_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_cfc12_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_sf6_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_oxy(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_tem(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_sal(nlons,nlats,ndeps,nmonths,nyears)
  ! variables to be calculated
  real :: mod_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_age(nlons,nlats,ndeps,nmonths,nyears)
  real :: obs_wma(nlons,nlats,ndeps,nmonths,nyears)
  real, dimension(:), allocatable :: amod_wma, amod_age, amod_oxy
  real, dimension(:), allocatable :: amod_ao2, amod_tem, amod_sal
  real, dimension(:), allocatable :: amod_rem, aobs_age, aobs_wma
  real, dimension(:), allocatable :: aobs_oxy, aobs_tem, aobs_sal
  real, dimension(:), allocatable :: ayear, amonth, adep, alat, alon
  ! temporary variables required during calculation
  ! for saving netcdf files with age data
  integer :: yeardim_id, monthdim_id, depdim_id, latdim_id, londim_id
  integer :: yearvar_id, monthvar_id, depvar_id, latvar_id, lonvar_id
  integer :: modwmavar_id, modagevar_id, modoxyvar_id, modao2var_id 
  integer :: modtemvar_id, modsalvar_id, modremvar_id, obsagevar_id
  integer :: obswmavar_id, obsoxyvar_id, obstemvar_id, obssalvar_id
  integer :: cntdim_id, dimids(5)
  

  !!! GET NETCDF DATA !!!

  print*, " *** "
  print*, " Reading model Water mass ages "
  call check( nf90_open(fmod_wma, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "Age", varid) )
  call check( nf90_get_var(ncid, varid, mod_wma(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model CFC-11 ages "
  call check( nf90_open(fmod_cfc11, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "CFC11age", varid) )
  call check( nf90_get_var(ncid, varid, mod_cfc11_age(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model CFC-12 ages "
  call check( nf90_open(fmod_cfc12, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "CFC12age", varid) )
  call check( nf90_get_var(ncid, varid, mod_cfc12_age(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model SF6 ages "
  call check( nf90_open(fmod_sf6, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "SF6age", varid) )
  call check( nf90_get_var(ncid, varid, mod_sf6_age(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model oxygen "
  call check( nf90_open(fmod_oxy, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "O2", varid) )
  call check( nf90_get_var(ncid, varid, mod_oxy(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check( nf90_inq_varid(ncid, "abioO2", varid) )
  call check( nf90_get_var(ncid, varid, mod_ao2(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model temperature "
  call check( nf90_open(fmod_tem, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "thetao", varid) )
  call check( nf90_get_var(ncid, varid, mod_tem(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model salinity "
  call check( nf90_open(fmod_sal, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "so", varid) )
  call check( nf90_get_var(ncid, varid, mod_sal(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading model remineralisation "
  call check( nf90_open(fmod_rem, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "REMIN", varid) )
  call check( nf90_get_var(ncid, varid, mod_rem(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check( nf90_inq_varid(ncid, "NITR", varid) )
  call check( nf90_get_var(ncid, varid, mod_nit(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  mod_rem(:,:,:,:,:) = (mod_rem(:,:,:,:,:)*(133/122.0) +     &
                        mod_nit(:,:,:,:,:)*2.0 )*365*86400*1e3
  
  print*, " *** "
  print*, " Reading observed CFC-11 ages "
  call check( nf90_open(fobs_cfc11, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "CFC11age", varid) )
  call check( nf90_get_var(ncid, varid, obs_cfc11_age(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading observed CFC-12 ages "
  call check( nf90_open(fobs_cfc12, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "CFC12age", varid) )
  call check( nf90_get_var(ncid, varid, obs_cfc12_age(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading observed SF6 ages "
  call check( nf90_open(fobs_sf6, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "SF6age", varid) )
  call check( nf90_get_var(ncid, varid, obs_sf6_age(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading observed oxygen "
  call check( nf90_open(fobs_oxy, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "oxygen", varid) )
  call check( nf90_get_var(ncid, varid, obs_oxy(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading observed temperature "
  call check( nf90_open(fobs_tem, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "temperature", varid) )
  call check( nf90_get_var(ncid, varid, obs_tem(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  print*, " *** "
  print*, " Reading observed salinity "
  call check( nf90_open(fobs_sal, nf90_nowrite, ncid) )
  call check( nf90_inq_varid(ncid, "salinity", varid) )
  call check( nf90_get_var(ncid, varid, obs_sal(:,:,:,:,:), &
              start=(/ 1,1,1,1,1 /), &
              count=(/ nlons, nlats, ndeps, nmonths, nyears /) ) )
  call check (nf90_close(ncid) )
  
  
  ! define years
  do m = 1,nyears
   years(m) = 1972.5 + (m-1)
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


  cnt = 0
  print*, " *** "
  print*, " Combine CFC/SF6 ages into one dataset "
  !!! Conduct calculation !!!
  do i = 1,nlons
  print*, " i ", i
   do j = 1,nlats
    do k = 1,ndeps
     do l = 1,nmonths
      do m = 1,nyears

        ! select ages using combination of CFCs and SF6
        if (years(m).lt.1990) then
          if (obs_cfc12_age(i,j,k,l,m).gt.0.0) then
            obs_age(i,j,k,l,m) = obs_cfc12_age(i,j,k,l,m)
            mod_age(i,j,k,l,m) = mod_cfc12_age(i,j,k,l,m)
            cnt = cnt + 1
          elseif (obs_cfc11_age(i,j,k,l,m).gt.0.0) then
            obs_age(i,j,k,l,m) = obs_cfc11_age(i,j,k,l,m)
            mod_age(i,j,k,l,m) = mod_cfc11_age(i,j,k,l,m)
            cnt = cnt + 1
          else
            obs_age(i,j,k,l,m) = -9999.0
            mod_age(i,j,k,l,m) = -9999.0
          endif
        else  ! years >= 1990
          if (obs_sf6_age(i,j,k,l,m).gt.0.0) then
            obs_age(i,j,k,l,m) = obs_sf6_age(i,j,k,l,m)
            mod_age(i,j,k,l,m) = mod_sf6_age(i,j,k,l,m)
            cnt = cnt + 1
          else
            obs_age(i,j,k,l,m) = -9999.0
            mod_age(i,j,k,l,m) = -9999.0
          endif
        endif

        ! if age is invalid, then remove data
        if (obs_age(i,j,k,l,m).le.0.0) then
          mod_wma(i,j,k,l,m) = -9999.0
          mod_oxy(i,j,k,l,m) = -9999.0
          mod_ao2(i,j,k,l,m) = -9999.0
          mod_tem(i,j,k,l,m) = -9999.0
          mod_sal(i,j,k,l,m) = -9999.0
          mod_rem(i,j,k,l,m) = -9999.0
          obs_oxy(i,j,k,l,m) = -9999.0
          obs_tem(i,j,k,l,m) = -9999.0
          obs_sal(i,j,k,l,m) = -9999.0
          obs_wma(i,j,k,l,m) = -9999.0
        else
          obs_wma(i,j,k,l,m) = (mod_wma(i,j,k,l,m)/mod_age(i,j,k,l,m)) * obs_age(i,j,k,l,m)
        endif

      enddo  ! m
     enddo  ! l
    enddo  ! k
   enddo  ! j
  enddo  ! i

  print*, "Number of valid age estimates = ", cnt

  allocate(alon(cnt))
  allocate(alat(cnt))
  allocate(adep(cnt))
  allocate(amonth(cnt))
  allocate(ayear(cnt))
  allocate(amod_wma(cnt))
  allocate(amod_age(cnt))
  allocate(amod_oxy(cnt))
  allocate(amod_ao2(cnt))
  allocate(amod_tem(cnt))
  allocate(amod_sal(cnt))
  allocate(amod_rem(cnt))
  allocate(aobs_age(cnt))
  allocate(aobs_wma(cnt))
  allocate(aobs_oxy(cnt))
  allocate(aobs_tem(cnt))
  allocate(aobs_sal(cnt))

  cnt2 = 0
  print*, " *** "
  print*, " Make 1D arrays for data saving "
  !!! Conduct calculation !!!
  do i = 1,nlons
  print*, " i ", i
   do j = 1,nlats
    do k = 1,ndeps
     do l = 1,nmonths
      do m = 1,nyears

        ! if age is invalid, then remove data
        if (obs_age(i,j,k,l,m).gt.0.0) then
          cnt2 = cnt2+1
          ! arrays for saving data as a table
          alon(cnt2) = lons(i)
          alat(cnt2) = lats(j)
          adep(cnt2) = deps(k)
          amonth(cnt2) = months(l)
          ayear(cnt2) = years(m)
          amod_wma(cnt2) = mod_wma(i,j,k,l,m)  
          amod_age(cnt2) = mod_age(i,j,k,l,m)  
          amod_oxy(cnt2) = mod_oxy(i,j,k,l,m)  
          amod_ao2(cnt2) = mod_ao2(i,j,k,l,m)  
          amod_tem(cnt2) = mod_tem(i,j,k,l,m)  
          amod_sal(cnt2) = mod_sal(i,j,k,l,m)  
          amod_rem(cnt2) = mod_rem(i,j,k,l,m)  
          aobs_age(cnt2) = obs_age(i,j,k,l,m)  
          aobs_wma(cnt2) = obs_wma(i,j,k,l,m)  
          aobs_oxy(cnt2) = obs_oxy(i,j,k,l,m)  
          aobs_tem(cnt2) = obs_tem(i,j,k,l,m)  
          aobs_sal(cnt2) = obs_sal(i,j,k,l,m)  
        endif

      enddo  ! m
     enddo  ! l
    enddo  ! k
   enddo  ! j
  enddo  ! i
  
  print*, " *** "
  print*, " Saving 1D array "
  call check( nf90_create("master.nc", nf90_clobber, ncid) )
  call check( nf90_def_dim(ncid, 'record', cnt, cntdim_id) )
  call check( nf90_def_var(ncid, 'lon', nf90_real, cntdim_id, lonvar_id) )
  call check( nf90_def_var(ncid, 'lat', nf90_real, cntdim_id, latvar_id) )
  call check( nf90_def_var(ncid, 'dep', nf90_real, cntdim_id, depvar_id) )
  call check( nf90_def_var(ncid, 'month', nf90_real, cntdim_id, monthvar_id) )
  call check( nf90_def_var(ncid, 'year', nf90_real, cntdim_id, yearvar_id) )
  call check( nf90_def_var(ncid, 'mod_wma', nf90_real, cntdim_id, modwmavar_id) )
  call check( nf90_def_var(ncid, 'mod_age', nf90_real, cntdim_id, modagevar_id) )
  call check( nf90_def_var(ncid, 'mod_oxy', nf90_real, cntdim_id, modoxyvar_id) )
  call check( nf90_def_var(ncid, 'mod_ao2', nf90_real, cntdim_id, modao2var_id) )
  call check( nf90_def_var(ncid, 'mod_tem', nf90_real, cntdim_id, modtemvar_id) )
  call check( nf90_def_var(ncid, 'mod_sal', nf90_real, cntdim_id, modsalvar_id) )
  call check( nf90_def_var(ncid, 'mod_rem', nf90_real, cntdim_id, modremvar_id) )
  call check( nf90_def_var(ncid, 'obs_age', nf90_real, cntdim_id, obsagevar_id) )
  call check( nf90_def_var(ncid, 'obs_wma', nf90_real, cntdim_id, obswmavar_id) )
  call check( nf90_def_var(ncid, 'obs_oxy', nf90_real, cntdim_id, obsoxyvar_id) )
  call check( nf90_def_var(ncid, 'obs_tem', nf90_real, cntdim_id, obstemvar_id) )
  call check( nf90_def_var(ncid, 'obs_sal', nf90_real, cntdim_id, obssalvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modwmavar_id, 'units', 'years') )
  call check( nf90_put_att(ncid, modagevar_id, 'units', 'years') )
  call check( nf90_put_att(ncid, modoxyvar_id, 'units', 'mmol / m3') )
  call check( nf90_put_att(ncid, modao2var_id, 'units', 'mmol / m3') )
  call check( nf90_put_att(ncid, modtemvar_id, 'units', 'degrees celcius') )
  call check( nf90_put_att(ncid, modsalvar_id, 'units', 'psu') )
  call check( nf90_put_att(ncid, modremvar_id, 'units', 'mol O2 / m3 / s') )
  call check( nf90_put_att(ncid, obsagevar_id, 'units', 'years') )
  call check( nf90_put_att(ncid, obswmavar_id, 'units', 'years') )
  call check( nf90_put_att(ncid, obsoxyvar_id, 'units', 'mmol / m3') )
  call check( nf90_put_att(ncid, obstemvar_id, 'units', 'degrees celcius') )
  call check( nf90_put_att(ncid, obssalvar_id, 'units', 'psu') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, ayear) )
  call check( nf90_put_var(ncid, monthvar_id, amonth) )
  call check( nf90_put_var(ncid, depvar_id, adep) )
  call check( nf90_put_var(ncid, latvar_id, alat) )
  call check( nf90_put_var(ncid, lonvar_id, alon) )
  call check( nf90_put_var(ncid, modwmavar_id, amod_wma) )
  call check( nf90_put_var(ncid, modagevar_id, amod_age) )
  call check( nf90_put_var(ncid, modoxyvar_id, amod_oxy) )
  call check( nf90_put_var(ncid, modao2var_id, amod_ao2) )
  call check( nf90_put_var(ncid, modtemvar_id, amod_tem) )
  call check( nf90_put_var(ncid, modsalvar_id, amod_sal) )
  call check( nf90_put_var(ncid, modremvar_id, amod_rem) )
  call check( nf90_put_var(ncid, obsagevar_id, aobs_age) )
  call check( nf90_put_var(ncid, obswmavar_id, aobs_wma) )
  call check( nf90_put_var(ncid, obsoxyvar_id, aobs_oxy) )
  call check( nf90_put_var(ncid, obstemvar_id, aobs_tem) )
  call check( nf90_put_var(ncid, obssalvar_id, aobs_sal) )
  call check( nf90_close(ncid) )

  
  print*, " *** "
  print*, " Saving model real water mass age values "
  call check( nf90_create("mod_wma.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_wma', nf90_real, dimids, modwmavar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modwmavar_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modwmavar_id, mod_wma, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving model age values "
  call check( nf90_create("mod_age.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_age', nf90_real, dimids, modagevar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modagevar_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modagevar_id, mod_age, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving model oxygen values "
  call check( nf90_create("mod_oxy.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_oxy', nf90_real, dimids, modoxyvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modoxyvar_id, 'units', 'mmol/m3') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modoxyvar_id, mod_oxy, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving model preformed oxygen values "
  call check( nf90_create("mod_ao2.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_ao2', nf90_real, dimids, modao2var_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modao2var_id, 'units', 'mmol/m3') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modao2var_id, mod_ao2, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving model temperature values "
  call check( nf90_create("mod_tem.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_tem', nf90_real, dimids, modtemvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modtemvar_id, 'units', 'degress celcius') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modtemvar_id, mod_tem, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving model salinity values "
  call check( nf90_create("mod_sal.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_sal', nf90_real, dimids, modsalvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modsalvar_id, 'units', 'psu') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modsalvar_id, mod_sal, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving model remineralisation values "
  call check( nf90_create("mod_rem.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'mod_rem', nf90_real, dimids, modremvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, modremvar_id, 'units', 'mol O2 / m3 / s') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, modremvar_id, mod_rem, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving observed age values "
  call check( nf90_create("obs_age.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'obs_age', nf90_real, dimids, obsagevar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, obsagevar_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, obsagevar_id, obs_age, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving observed wma values "
  call check( nf90_create("obs_wma.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'obs_wma', nf90_real, dimids, obswmavar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, obswmavar_id, 'units', 'years') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, obswmavar_id, obs_wma, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving observed oxygen values "
  call check( nf90_create("obs_oxy.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'obs_oxy', nf90_real, dimids, obsoxyvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, obsoxyvar_id, 'units', 'mmol/m3') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, obsoxyvar_id, obs_oxy, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving observed temperature values "
  call check( nf90_create("obs_tem.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'obs_tem', nf90_real, dimids, obstemvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, obstemvar_id, 'units', 'degress celcius') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, obstemvar_id, obs_tem, &
              start=(/ 1, 1, 1, 1, 1/), count=(/ nlons, nlats, ndeps, nmonths, nyears /)) )
  call check( nf90_close(ncid) )

  print*, " *** "
  print*, " Saving observed salinity values "
  call check( nf90_create("obs_sal.nc", nf90_clobber, ncid) )
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
  call check( nf90_def_var(ncid, 'obs_sal', nf90_real, dimids, obssalvar_id) )
  call check( nf90_put_att(ncid, yearvar_id, 'units', 'year (common era)') )
  call check( nf90_put_att(ncid, monthvar_id, 'units', 'months') )
  call check( nf90_put_att(ncid, depvar_id, 'units', 'metres') )
  call check( nf90_put_att(ncid, latvar_id, 'units', 'degrees_north') )
  call check( nf90_put_att(ncid, lonvar_id, 'units', 'degrees_east') )
  call check( nf90_put_att(ncid, obssalvar_id, 'units', 'psu') )
  call check( nf90_enddef(ncid) )
  call check( nf90_put_var(ncid, yearvar_id, years) )
  call check( nf90_put_var(ncid, monthvar_id, months) )
  call check( nf90_put_var(ncid, depvar_id, deps) )
  call check( nf90_put_var(ncid, latvar_id, lats) )
  call check( nf90_put_var(ncid, lonvar_id, lons) )
  call check( nf90_put_var(ncid, obssalvar_id, obs_sal, &
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
