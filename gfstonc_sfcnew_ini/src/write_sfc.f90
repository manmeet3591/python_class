! theia:
! gfortran -O3 -march=native -O3 -fPIC -c kinds.f90
! gfortran -O3 -march=native -O3 -fPIC -c sfcio_module.f90
! f2py -c write_sfc.f90 -m write_sfc --fcompiler=gnu95 kinds.o sfcio_module.o
!subroutine write_header(filename, sfchead)
!  use sfcio_module, only: sfcio_sclose,sfcio_head,sfcio_swhead,&
!  sfcio_swohdc,sfcio_aldata,sfcio_data,sfcio_swopen,sfcio_swdata,sfcio_axdata
!  implicit none
!  integer :: nlons,nlats,idate(4),lsoil
!  real :: fhour
!  character(len=500), intent(in) :: filename
!  type(sfcio_head), intent(in) :: sfchead
!  integer lu,iret
!  lu = 7
!  call sfcio_swopen(lu,trim(filename),iret)
!  if (iret .ne. 0) then
!     print *,'error opening ',trim(filename),iret
!     stop
!  endif
!  call sfcio_swhead(lu,sfchead,iret)
!  if (iret .ne. 0) then
!     print *,'error reading header from ',trim(filename),iret
!     stop
!  else
!     nlons = sfchead%lonb
!     nlats = sfchead%latb
!     idate = sfchead%idate
!     lsoil = sfchead%lsoil
!     fhour = sfchead%fhour
!  endif 
!  call sfcio_sclose(lu,iret)
!end subroutine write_header


!subroutine read_header(filename, nlons, nlats, lsoil, idate, fhour)
!  use sfcio_module, only: sfcio_sclose,sfcio_head,sfcio_srhead,&
!  sfcio_srohdc,sfcio_aldata,sfcio_data,sfcio_sropen,sfcio_srdata,sfcio_axdata
!  implicit none
!  integer, intent(out) :: nlons,nlats,idate(4),lsoil
!  real, intent(out) :: fhour
!  character(len=500), intent(in) :: filename
!  type(sfcio_head) sfchead
!  integer lu,iret
!  lu = 7
!  call sfcio_sropen(lu,trim(filename),iret)
!  if (iret .ne. 0) then
!     print *,'error opening ',trim(filename),iret
!     stop
!  endif
!  call sfcio_srhead(lu,sfchead,iret)
!  if (iret .ne. 0) then
!     print *,'error reading header from ',trim(filename),iret
!     stop
!  else
!     nlons = sfchead%lonb
!     nlats = sfchead%latb
!     idate = sfchead%idate
!     lsoil = sfchead%lsoil
!     fhour = sfchead%fhour
!  endif 
!  call sfcio_sclose(lu,iret)
!end subroutine read_header

subroutine write_griddata(filename, nlons, nlats, lsoil, idate, fhour,grids2d,grids3d)
  use sfcio_module, only: sfcio_sclose,sfcio_swohdc,sfcio_head,sfcio_swhead,&
  sfcio_srohdc,sfcio_aldata,sfcio_data,sfcio_swopen,sfcio_swdata,sfcio_axdata
  use kinds, only: r_kind
  implicit none
  integer, parameter :: n2d = 32
  integer, parameter :: n3d = 3
  integer, intent(in) :: nlons,nlats,lsoil, idate(4)
  real(r_kind), intent(in),dimension(nlons,nlats,n2d) :: grids2d
  real(r_kind), intent(in), dimension(nlons,nlats,lsoil,n3d) :: grids3d
  real, intent(in) :: fhour
  character(len=*), intent(in) :: filename
  integer lu,iret
  type(sfcio_head) :: sfchead
  type(sfcio_data) :: sfcdata
  integer,parameter:: sfcio_intkind=4,sfcio_realkind=4
  integer,parameter:: sfcio_lhead1=32
 ! integer(sfcio_intkind),allocatable:: lpl(:)
  real(sfcio_realkind):: zsoil(4)
  character(sfcio_lhead1) :: clabsfc='                                '

!  allocate(lpl(int(nlats/2)))


  lu = 8
  !print *,"Entered write_griddata", filename, nlons
  sfchead%lonb      =    nlons
  print *,"test1 lonb filename lsoil", sfchead%lonb, filename, lsoil
  sfchead%latb      =    nlats
  sfchead%idate     =    idate
  sfchead%lsoil     =    lsoil
  sfchead%fhour     =    fhour
  sfchead%ivs       =    200501
  if (nlats.eq.94) then
  sfchead%lpl       =    (/30,  30,  30,  40,  48,  56,  60,  72,  72,  80,  90,  90,       &
						   96, 110, 110, 120, 120, 128, 144, 144, 144, 144, 154, 160,       &
						  160, 168, 168, 180, 180, 180, 180, 180, 180, 192, 192, 192,       &
						  192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192/)
  else if (nlats.eq.190) then
  sfchead%lpl       =    (/30,   30,   36,   48,   56,   60,   72,   72,   80,   90,      &
						   96,  110,  110,  120,  120,  128,  144,  144,  154,  160,      &
						   160,  180,  180,  180,  192,  192,  210,  210,  220,  220,      &
						   240,  240,  240,  240,  240,  252,  256,  280,  280,  280,      &
						   280,  288,  288,  288,  288,  308,  308,  308,  320,  320,      &
						   320,  320,  330,  330,  360,  360,  360,  360,  360,  360,      &
						   360,  360,  360,  360,  360,  360,  384,  384,  384,  384,      &
						   384,  384,  384,  384,  384,  384,  384,  384,  384,  384,      &
						   384,  384,  384,  384,  384,  384,  384,  384,  384,  384,      &
						   384,  384,  384,  384,  384/)

  end if

  zsoil     =   (/-0.1,-0.4,-1.0,-2.0/) 
  sfchead%zsoil     = zsoil 
  sfchead%clabsfc   = clabsfc

  !print *,"test3 fhour shape(grids2d)", sfchead%fhour, shape(grids2d)
  allocate(sfcdata%slmsk(nlons, nlats))
  sfcdata%slmsk   =  grids2d(:,:,1)   
  allocate(sfcdata%orog(nlons, nlats))
  sfcdata%orog    =  grids2d(:,:,2)   
  allocate(sfcdata%tsea(nlons, nlats))
  sfcdata%tsea    =  grids2d(:,:,3)   
  allocate(sfcdata%sheleg(nlons, nlats))
  sfcdata%sheleg  =  grids2d(:,:,4)   
  allocate(sfcdata%tg3(nlons, nlats))
  sfcdata%tg3     =  grids2d(:,:,5)   
  allocate(sfcdata%zorl(nlons, nlats))
  sfcdata%zorl    =  grids2d(:,:,6)   
  allocate(sfcdata%alvsf(nlons, nlats))
  sfcdata%alvsf   =  grids2d(:,:,7)   
  allocate(sfcdata%alvwf(nlons, nlats))
  sfcdata%alvwf   =  grids2d(:,:,8)   
  allocate(sfcdata%alnsf(nlons, nlats))
  sfcdata%alnsf   =  grids2d(:,:,9)   
  allocate(sfcdata%alnwf(nlons, nlats))
  sfcdata%alnwf   =  grids2d(:,:,10)  
  allocate(sfcdata%vfrac(nlons, nlats))
  sfcdata%vfrac   =  grids2d(:,:,11)  
  allocate(sfcdata%canopy(nlons, nlats))
  sfcdata%canopy  =  grids2d(:,:,12)  
  allocate(sfcdata%f10m(nlons, nlats))
  sfcdata%f10m    =  grids2d(:,:,13)  
  allocate(sfcdata%t2m(nlons, nlats))
  sfcdata%t2m     =  grids2d(:,:,14)  
  allocate(sfcdata%q2m(nlons, nlats))
  sfcdata%q2m     =  grids2d(:,:,15)  
  allocate(sfcdata%vtype(nlons, nlats))
  sfcdata%vtype   =  grids2d(:,:,16)  
  allocate(sfcdata%stype(nlons, nlats))
  sfcdata%stype   =  grids2d(:,:,17)  
  allocate(sfcdata%facsf(nlons, nlats))
  sfcdata%facsf   =  grids2d(:,:,18)  
  allocate(sfcdata%facwf(nlons, nlats))
  sfcdata%facwf   =  grids2d(:,:,19)  
  allocate(sfcdata%uustar(nlons, nlats))
  sfcdata%uustar  =  grids2d(:,:,20)  
  allocate(sfcdata%ffmm(nlons, nlats))
  sfcdata%ffmm    =  grids2d(:,:,21)  
  allocate(sfcdata%ffhh(nlons, nlats))
  sfcdata%ffhh    =  grids2d(:,:,22)  
  allocate(sfcdata%hice(nlons, nlats))
  sfcdata%hice    =  grids2d(:,:,23)  
  allocate(sfcdata%fice(nlons, nlats))
  sfcdata%fice    =  grids2d(:,:,24)  
  allocate(sfcdata%tisfc(nlons, nlats))
  sfcdata%tisfc   =  grids2d(:,:,25)  
  allocate(sfcdata%tprcp(nlons, nlats))
  sfcdata%tprcp   =  grids2d(:,:,26)  
  allocate(sfcdata%srflag(nlons, nlats))
  sfcdata%srflag  =  grids2d(:,:,27)  
  allocate(sfcdata%snwdph(nlons, nlats))
  sfcdata%snwdph  =  grids2d(:,:,28)  
  allocate(sfcdata%shdmin(nlons, nlats))
  sfcdata%shdmin  =  grids2d(:,:,29)  
  allocate(sfcdata%shdmax(nlons, nlats))
  sfcdata%shdmax  =  grids2d(:,:,30)  
  allocate(sfcdata%slope(nlons, nlats))
  sfcdata%slope   =  grids2d(:,:,31)  
  allocate(sfcdata%snoalb(nlons, nlats))
  sfcdata%snoalb  =  grids2d(:,:,32)  
  allocate(sfcdata%stc(nlons, nlats, lsoil))
  sfcdata%stc     =  grids3d(:,:,:,1) 
  allocate(sfcdata%smc(nlons, nlats, lsoil))
  sfcdata%smc     =  grids3d(:,:,:,2) 
  allocate(sfcdata%slc(nlons, nlats, lsoil))
  sfcdata%slc     =  grids3d(:,:,:,3) 
 print *,"test 2 sfcdata%smc(1,1,1) ",sfcdata%smc(1,1,1)

  call sfcio_swohdc(lu,trim(filename),sfchead,sfcdata,iret)


  if (iret .ne. 0) then
    print *,'error writing ',trim(filename),iret
    stop
  endif


end subroutine write_griddata


!subroutine strtoarr(strin, chararr, n_str)
!  integer, intent(in) :: n_str
!  character(len=n_str), intent(in) :: strin
!  integer, intent(out) ::  chararr(n_str)
!  chararr = 32
!  do j=1,len_trim(trim(adjustl(strin)))
!     chararr(j) = ichar(strin(j:j))
!  enddo
!end subroutine strtoarr
