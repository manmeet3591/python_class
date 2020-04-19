! theia:
! gfortran -O3 -march=native -O3 -fPIC -c kinds.f90
! gfortran -O3 -march=native -O3 -fPIC -c sfcio_module.f90
! f2py -c read_sfc.f90 -m read_sfc --fcompiler=gnu95 kinds.o sfcio_module.o
subroutine read_header(filename, nlons, nlats, lsoil, idate, fhour)
  use sfcio_module, only: sfcio_sclose,sfcio_head,sfcio_srhead,&
  sfcio_srohdc,sfcio_aldata,sfcio_data,sfcio_sropen,sfcio_srdata,sfcio_axdata
  implicit none
  integer, intent(out) :: nlons,nlats,idate(4),lsoil
  real, intent(out) :: fhour
  character(len=500), intent(in) :: filename
  type(sfcio_head) sfchead
  integer lu,iret
  lu = 7
  call sfcio_sropen(lu,trim(filename),iret)
  if (iret .ne. 0) then
     print *,'error opening ',trim(filename),iret
     stop
  endif
  call sfcio_srhead(lu,sfchead,iret)
  if (iret .ne. 0) then
     print *,'error reading header from ',trim(filename),iret
     stop
  else
     nlons = sfchead%lonb
     nlats = sfchead%latb
     idate = sfchead%idate
     lsoil = sfchead%lsoil
     fhour = sfchead%fhour
  endif 
  call sfcio_sclose(lu,iret)
end subroutine read_header

subroutine read_griddata(filename, nlons, nlats, lsoil,&
                grids2d,grids2d_desc,grids2d_name,grids3d,grids3d_desc,grids3d_name)
  use sfcio_module, only: sfcio_sclose,sfcio_swohdc,sfcio_head,sfcio_srhead,&
  sfcio_srohdc,sfcio_aldata,sfcio_data,sfcio_sropen,sfcio_srdata,sfcio_axdata
  use kinds, only: r_kind
  implicit none
  integer, parameter :: n_str1 = 72
  integer, parameter :: n_str2 = 8
  integer, parameter :: n2d = 32
  integer, parameter :: n3d = 3
  character(len=n_str1) :: desc
  character(len=n_str2) :: name
  real(r_kind), intent(out),dimension(nlons,nlats,n2d) :: grids2d
  integer, intent(out), dimension(n2d,n_str1) :: grids2d_desc
  integer, intent(out), dimension(n2d,n_str2) :: grids2d_name
  real(r_kind), intent(out), dimension(nlons,nlats,lsoil,n3d) :: grids3d
  integer, intent(out), dimension(n3d,n_str1) :: grids3d_desc
  integer, intent(out), dimension(n3d,n_str1) :: grids3d_name
  integer, intent(in) :: nlons,nlats,lsoil
  character(len=500), intent(in) :: filename
  integer lu,iret
  type(sfcio_head) :: sfchead
  type(sfcio_data) :: sfcdata
  lu = 7
  !print *,trim(filename)
  call sfcio_srohdc(lu,trim(filename),sfchead,sfcdata,iret)
  if (iret .ne. 0) then
    print *,'error reading ',trim(filename),iret
    stop
  endif
  desc = 'sea-land-ice mask (0-sea, 1-land, 2-ice)'
  call strtoarr(desc, grids2d_desc(1,:), n_str1)
  name = 'slmsk'
  call strtoarr(name, grids2d_name(1,:), n_str2)
  grids2d(:,:,1) = sfcdata%slmsk
  desc = 'surface orography in m'
  call strtoarr(desc, grids2d_desc(2,:), n_str1)
  name = 'orog'
  call strtoarr(name, grids2d_name(2,:), n_str2)
  grids2d(:,:,2) = sfcdata%orog
  desc = 'sst in K'
  call strtoarr(desc, grids2d_desc(3,:), n_str1)
  name = 'tsea'
  call strtoarr(name, grids2d_name(3,:), n_str2)
  grids2d(:,:,3) = sfcdata%tsea
  desc = 'snow depth in m'
  call strtoarr(desc, grids2d_desc(4,:), n_str1)
  name = 'sheleg'
  call strtoarr(name, grids2d_name(4,:), n_str2)
  grids2d(:,:,4) = sfcdata%sheleg
  desc = 'deep soil temp in K'
  call strtoarr(desc, grids2d_desc(5,:), n_str1)
  name = 'tg3'
  call strtoarr(name, grids2d_name(5,:), n_str2)
  grids2d(:,:,5) = sfcdata%tg3
  desc = 'roughness length in cm'
  call strtoarr(desc, grids2d_desc(6,:), n_str1)
  name = 'zorl'
  call strtoarr(name, grids2d_name(6,:), n_str2)
  grids2d(:,:,6) = sfcdata%zorl
  desc = 'albedo for visible scattered (fraction)'
  call strtoarr(desc, grids2d_desc(7,:), n_str1)
  name = 'alvsf'
  call strtoarr(name, grids2d_name(7,:), n_str2)
  grids2d(:,:,7) = sfcdata%alvsf
  desc = 'albedo for visible beam (fraction)'
  call strtoarr(desc, grids2d_desc(8,:), n_str1)
  name = 'alvwf'
  call strtoarr(name, grids2d_name(8,:), n_str2)
  grids2d(:,:,8) = sfcdata%alvwf
  desc = 'albedo for near-IR scattered (fraction)'
  call strtoarr(desc, grids2d_desc(9,:), n_str1)
  name = 'alnsf'
  call strtoarr(name, grids2d_name(9,:), n_str2)
  grids2d(:,:,9) = sfcdata%alnsf
  desc = 'albedo for near-IR beam (fraction)'
  call strtoarr(desc, grids2d_desc(10,:), n_str1)
  name = 'alnwf'
  call strtoarr(name, grids2d_name(10,:), n_str2)
  grids2d(:,:,10) = sfcdata%alnwf
  desc = 'vegetation fraction (fraction)'
  call strtoarr(desc, grids2d_desc(11,:), n_str1)
  name = 'vfrac'
  call strtoarr(name, grids2d_name(11,:), n_str2)
  grids2d(:,:,11) = sfcdata%vfrac
  desc = 'canopy water in m'
  call strtoarr(desc, grids2d_desc(12,:), n_str1)
  name = 'canopy'
  call strtoarr(name, grids2d_name(12,:), n_str2)
  grids2d(:,:,12) = sfcdata%canopy
  desc = '10-meter wind speed over lowest model wind speed'
  call strtoarr(desc, grids2d_desc(13,:), n_str1)
  name = 'f10m'
  call strtoarr(name, grids2d_name(13,:), n_str2)
  grids2d(:,:,13) = sfcdata%f10m
  desc = '2-meter temp in K'
  call strtoarr(desc, grids2d_desc(14,:), n_str1)
  name = 't2m'
  call strtoarr(name, grids2d_name(14,:), n_str2)
  grids2d(:,:,14) = sfcdata%t2m
  desc = '2-meter specific humidity (kg/kg)'
  call strtoarr(desc, grids2d_desc(15,:), n_str1)
  name = 'q2m'
  call strtoarr(name, grids2d_name(15,:), n_str2)
  grids2d(:,:,15) = sfcdata%q2m
  desc = 'vegetation type in integer 1-13'
  call strtoarr(desc, grids2d_desc(16,:), n_str1)
  name = 'vtype'
  call strtoarr(name, grids2d_name(16,:), n_str2)
  grids2d(:,:,16) = sfcdata%vtype
  desc = 'soil type in integer 1-9'
  call strtoarr(desc, grids2d_desc(17,:), n_str1)
  name = 'stype'
  call strtoarr(name, grids2d_name(17,:), n_str2)
  grids2d(:,:,17) = sfcdata%stype
  desc = 'fractional coverage with strong cosz dependency'
  call strtoarr(desc, grids2d_desc(18,:), n_str1)
  name = 'facsf'
  call strtoarr(name, grids2d_name(18,:), n_str2)
  grids2d(:,:,18) = sfcdata%facsf
  desc = 'fractional coverage with weak cosz dependency'  
  call strtoarr(desc, grids2d_desc(19,:), n_str1)
  name = 'facwf'
  call strtoarr(name, grids2d_name(19,:), n_str2)
  grids2d(:,:,19) = sfcdata%facwf
  desc = 'surface layer friction velocity (m/s)'
  call strtoarr(desc, grids2d_desc(20,:), n_str1)
  name = 'uustar'
  call strtoarr(name, grids2d_name(20,:), n_str2)
  grids2d(:,:,20) = sfcdata%uustar
  desc = 'F_m parameter from PBL scheme = LOG((Z0MAX(I)+Z1(I)) / Z0MAX(I))'
  call strtoarr(desc, grids2d_desc(21,:), n_str1)
  name = 'ffmm'
  call strtoarr(name, grids2d_name(21,:), n_str2)
  grids2d(:,:,21) = sfcdata%ffmm
  desc = 'F_h parameter from PBL scheme = LOG((ZTMAX(I)+Z1(I)) / ZTMAX(I))'
  call strtoarr(desc, grids2d_desc(22,:), n_str1)
  name = 'ffhh'
  call strtoarr(name, grids2d_name(22,:), n_str2)
  grids2d(:,:,22) = sfcdata%ffhh
  desc = 'sea-ice thickness'
  call strtoarr(desc, grids2d_desc(23,:), n_str1)
  name = 'hice'
  call strtoarr(name, grids2d_name(23,:), n_str2)
  grids2d(:,:,23) = sfcdata%hice
  desc = 'sea-ice concentration'
  call strtoarr(desc, grids2d_desc(24,:), n_str1)
  name = 'fice'
  call strtoarr(name, grids2d_name(24,:), n_str2)
  grids2d(:,:,24) = sfcdata%fice
  desc = 'sea-ice temperature'
  call strtoarr(desc, grids2d_desc(25,:), n_str1)
  name = 'tisfc'
  call strtoarr(name, grids2d_name(25,:), n_str2)
  grids2d(:,:,25) = sfcdata%tisfc
  desc = 'accumulated total precipitation (kg/m2) '
  call strtoarr(desc, grids2d_desc(26,:), n_str1)
  name = 'tprcp'
  call strtoarr(name, grids2d_name(26,:), n_str2)
  grids2d(:,:,26) = sfcdata%tprcp
  desc = 'snow/rain flag for precipitation (0 rain, 1 snow)'
  call strtoarr(desc, grids2d_desc(27,:), n_str1)
  name = 'srflag'
  call strtoarr(name, grids2d_name(27,:), n_str2)
  grids2d(:,:,27) = sfcdata%srflag
  desc = 'actual snow depth (mm) over land/sea ice'
  call strtoarr(desc, grids2d_desc(28,:), n_str1)
  name = 'snwdph'
  call strtoarr(name, grids2d_name(28,:), n_str2)
  grids2d(:,:,28) = sfcdata%snwdph
  desc = 'minimum areal coverage of green veg'
  call strtoarr(desc, grids2d_desc(29,:), n_str1)
  name = 'shdmin'
  call strtoarr(name, grids2d_name(29,:), n_str2)
  grids2d(:,:,29) = sfcdata%shdmin
  desc = 'maximum areal coverage of green veg'
  call strtoarr(desc, grids2d_desc(30,:), n_str1)
  name = 'shdmax'
  call strtoarr(name, grids2d_name(30,:), n_str2)
  grids2d(:,:,30) = sfcdata%shdmax
  desc = 'integer class of sfc slope '
  call strtoarr(desc, grids2d_desc(31,:), n_str1)
  name = 'slope'
  call strtoarr(name, grids2d_name(31,:), n_str2)
  grids2d(:,:,31) = sfcdata%slope
  desc = 'maximum (deep) snow albedo'
  call strtoarr(desc, grids2d_desc(32,:), n_str1)
  name = 'snoalb'
  call strtoarr(name, grids2d_name(32,:), n_str2)
  grids2d(:,:,32) = sfcdata%snoalb
  desc = 'soil temperature in K'
  call strtoarr(desc, grids3d_desc(1,:), n_str1)
  name = 'stc'
  call strtoarr(name, grids3d_name(1,:), n_str2)
  grids3d(:,:,:,1) = sfcdata%stc
  desc = 'soil volumetric water content (fraction)'
  call strtoarr(desc, grids3d_desc(2,:), n_str1)
  name = 'smc'
  call strtoarr(name, grids3d_name(2,:), n_str2)
  grids3d(:,:,:,2) = sfcdata%smc
  desc = 'liquid soil moisture content (fraction)'
  call strtoarr(desc, grids3d_desc(3,:), n_str1)
  name = 'slc'
  call strtoarr(name, grids3d_name(3,:), n_str2)
  grids3d(:,:,:,3) = sfcdata%slc
end subroutine read_griddata
subroutine strtoarr(strin, chararr, n_str)
  integer, intent(in) :: n_str
  character(len=n_str), intent(in) :: strin
  integer, intent(out) ::  chararr(n_str)
  chararr = 32
  do j=1,len_trim(trim(adjustl(strin)))
     chararr(j) = ichar(strin(j:j))
  enddo
end subroutine strtoarr
