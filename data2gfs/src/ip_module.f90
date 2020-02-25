module ip_module
  use type_module
  use constant_module, only: math_pi
  implicit none

  private

  public :: kgds2kgds,makekgds
  
  interface kgds2kgds
     module procedure kgds2kgds_scalar_2d, kgds2kgds_scalar_3d, &
          &           kgds2kgds_vector_2d, kgds2kgds_vector_3d
  end interface kgds2kgds
  
contains

  subroutine kgds2kgds_scalar_2d(ip,fldin,kgdsi,kgdso,fldout,maxwv)
    integer(kind=i4b),intent(in) :: ip   !(0=bilinear,1=bicubic,2=neighbor)
    real(kind=sp),intent(in)     :: fldin(:,:)
    integer(kind=i4b),intent(in) :: kgdsi(200),kgdso(200)
    integer(kind=i4b),intent(in),optional :: maxwv
    real(kind=sp),intent(inout)  :: fldout(:,:)

    integer(kind=i4b)  :: xi,yi,xo,yo,km,no,iret
    integer(kind=i4b)  :: ipopt(20),ibi(1),ibo(1)
    real(kind=i4b)     :: rlat(size(fldin)),rlon(size(fldout))
    logical            :: li(size(fldin)),lo(size(fldout))

    xi=size(fldin,1)
    yi=size(fldin,2)
    xo=size(fldout,1)
    yo=size(fldout,2)
    km=1    
    ipopt = -1   
    if (ip == 0) then
       write(*,*) 'bilinear interpolation'
    else if (ip == 1) then 
       write(*,*) 'bicubic interpolation'
    else if (ip == 4) then
       ipopt(1) = 0! for triangular
       if (present(maxwv) .and. maxwv>0) then
          ipopt(2) = maxwv
       else
          if (kgdsi(1)==0) then  ! latitude/longitude grid
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-3)/2) !linear grid
          else if (kgdsi(1)==4) then ! gausian grid
             ipopt(2) = min(kgdsi(2)/2-1, kgdsi(3)-1)     !linear grid
          else if (kgdsi(1)==256) then ! latitude/longitude grid excluding pole
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-1)/2) !linear grid
          else
             ipopt(2) = -1
          end if
       end if
       write(*,*) 'spectral interpolation. max wave number=',ipopt(2)
    else
       write(*,*) 'error! invalid iptypes. ip=', ip
       call abort
    end if
    ibi   = 0
    li    = .false.
    call ipolates(ip,ipopt,kgdsi,kgdso,xi*yi,xo*yo,km,ibi,li,fldin,  &
         &        no,rlat,rlon,ibo,lo,fldout,iret)

  end subroutine kgds2kgds_scalar_2d
  subroutine kgds2kgds_scalar_3d(ip,fldin,kgdsi,kgdso,fldout,maxwv)
    integer(kind=i4b),intent(in) :: ip   !(0=bilinear,1=bicubic,4=spectral)
    real(kind=sp),intent(in)     :: fldin(:,:,:)
    integer(kind=i4b),intent(in) :: kgdsi(200),kgdso(200)
    integer(kind=i4b),intent(in),optional :: maxwv
    real(kind=sp),intent(inout)  :: fldout(:,:,:)

    integer(kind=i4b)  :: xi,yi,xo,yo,km,no,iret
    integer(kind=i4b)  :: ipopt(20),ibi(size(fldin,3)),ibo(size(fldout,3))
    real(kind=i4b)     :: rlat(size(fldout,1)*size(fldout,2)),rlon(size(fldout,1)*size(fldout,2))
    logical            :: li(size(fldin)),lo(size(fldout))

    xi=size(fldin,1)
    yi=size(fldin,2)
    xo=size(fldout,1)
    yo=size(fldout,2)
    km=size(fldout,3)
 
    ipopt = -1   
    if (ip == 0) then
       write(*,*) 'bilinear interpolation'
    else if (ip == 1) then 
       write(*,*) 'bicubic interpolation'
    else if (ip == 4) then
       ipopt(1) = 0! for triangular
       if (present(maxwv) .and. maxwv>0) then
          ipopt(2) = maxwv
       else
          if (kgdsi(1)==0) then  ! latitude/longitude grid
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-3)/2) !linear grid
          else if (kgdsi(1)==4) then ! gausian grid
             ipopt(2) = min(kgdsi(2)/2-1, kgdsi(3)-1)     !linear grid
          else if (kgdsi(1)==256) then ! latitude/longitude grid excluding pole
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-1)/2) !linear grid
          else
             ipopt(2) = -1
          end if
       end if
       write(*,*) 'spectral interpolation. max wave number=',ipopt(2)
    else
       write(*,*) 'error! invalid iptypes. ip=', ip
       call abort
    end if
    ibi   = 0
    li    = .false.
    call ipolates(ip,ipopt,kgdsi,kgdso,xi*yi,xo*yo,km,ibi,li,fldin,  &
         &        no,rlat,rlon,ibo,lo,fldout,iret)

  end subroutine kgds2kgds_scalar_3d

  subroutine kgds2kgds_vector_2d(ip,uin,vin,kgdsi,kgdso,uout,vout,maxwv)
    integer(kind=i4b),intent(in) :: ip   !(0=bilinear,1=bicubic,2=neighbor)
    real(kind=sp),intent(in)     :: uin(:,:),vin(:,:)
    integer(kind=i4b),intent(in) :: kgdsi(200),kgdso(200)
    integer(kind=i4b),intent(in),optional :: maxwv
    real(kind=sp),intent(inout)  :: uout(:,:),vout(:,:)

    integer(kind=i4b)  :: xi,yi,xo,yo,km,no,iret
    integer(kind=i4b)  :: ipopt(20),ibi(1),ibo(1)
    real(kind=i4b)     :: rlat(size(uin)),rlon(size(uout))
    logical            :: li(size(uin)),lo(size(uout))

    xi=size(uin,1)
    yi=size(uin,2)
    xo=size(uout,1)
    yo=size(uout,2)
    km=1    
    ibi   = 0
    li    = .false.
    ipopt = -1   
    if (ip == 0) then
       write(*,*) 'bilinear interpolation'
    else if (ip == 1) then 
       write(*,*) 'bicubic interpolation'
    else if (ip == 4) then
       ipopt(1) = 0 ! for triangular
       if (present(maxwv) .and. maxwv>0) then
          ipopt(2) = maxwv
       else
          if (kgdsi(1)==0) then  ! latitude/longitude grid
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-3)/2) !linear grid
          else if (kgdsi(1)==4) then ! gausian grid
             ipopt(2) = min(kgdsi(2)/2-1, kgdsi(3)-1)     !linear grid
          else if (kgdsi(1)==256) then ! latitude/longitude grid excluding pole
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-1)/2) !linear grid
          else
             ipopt(2) = -1
          end if
       end if
       write(*,*) 'spectral interpolation. max wave number=',ipopt(2)
    else
       write(*,*) 'error! invalid iptypes. ip=', ip
       call abort
    end if
    call ipolatev(ip,ipopt,kgdsi,kgdso,xi*yi,xo*yo,km,ibi,li,uin,vin,  &
         &        no,rlat,rlon,ibo,lo,uout,vout,iret)

  end subroutine kgds2kgds_vector_2d
  subroutine kgds2kgds_vector_3d(ip,uin,vin,kgdsi,kgdso,uout,vout,maxwv)
    integer(kind=i4b),intent(in) :: ip   !(0=bilinear,1=bicubic,2=neighbor)
    real(kind=sp),intent(in)     :: uin(:,:,:),vin(:,:,:)
    integer(kind=i4b),intent(in) :: kgdsi(200),kgdso(200)
    integer(kind=i4b),intent(in),optional :: maxwv
    real(kind=sp),intent(inout)  :: uout(:,:,:),vout(:,:,:)

    integer(kind=i4b)  :: xi,yi,xo,yo,km,no,iret
    integer(kind=i4b)  :: ipopt(20),ibi(size(uin,3)),ibo(size(uout,3))
    real(kind=i4b)     :: rlat(size(uout,1)*size(uout,2)),rlon(size(uout,1)*size(uout,2)), &
         &                crot(size(uout,1)*size(uout,2)),srot(size(uout,1)*size(uout,2))
    logical            :: li(size(uin)),lo(size(uout))

    xi=size(uin,1)
    yi=size(uin,2)
    xo=size(uout,1)
    yo=size(uout,2)
    km=size(uout,3)
    ibi   = 0
    li    = .false.
    ipopt = -1   
    if (ip == 0) then
       write(*,*) 'bilinear interpolation'
    else if (ip == 1) then 
       write(*,*) 'bicubic interpolation'
    else if (ip == 4) then
       ipopt(1) = 0 ! for triangular
       if (present(maxwv) .and. maxwv>0) then
          ipopt(2) = maxwv
       else
          if (kgdsi(1)==0) then  ! latitude/longitude grid
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-3)/2) !linear grid
          else if (kgdsi(1)==4) then ! gausian grid
             ipopt(2) = min(kgdsi(2)/2-1, kgdsi(3)-1)     !linear grid
          else if (kgdsi(1)==256) then ! latitude/longitude grid excluding pole
             ipopt(2) = min(kgdsi(2)/2-1, (kgdsi(3)-1)/2) !linear grid
          else
             ipopt(2) = -1
          end if
       end if
       write(*,*) 'spectral interpolation. max wave number=',ipopt(2)
    else
       write(*,*) 'error! invalid iptypes. ip=', ip
       call abort
    end if
    call ipolatev(ip,ipopt,kgdsi,kgdso,xi*yi,xo*yo,km,ibi,li,uin,vin,  &
         &        no,rlat,rlon,crot,srot,ibo,lo,uout,vout,iret)
  end subroutine kgds2kgds_vector_3d

  subroutine makekgds(grid,in,jn,yrev,kgds)
    character(len=2), intent(in)   :: grid 
    integer(kind=i4b), intent(in)  :: in, jn
    logical, intent(in)            :: yrev
    !yrev=.true. latitude scans negaticely
    integer(kind=i4b),intent(out)  :: kgds(200)
    real(kind=sp) :: slat(jn),wlat(jn)
    real(kind=sp),parameter :: pi=math_pi

    kgds=0
    if (grid == 'll') then !latitude-londitude equally spaced grid
       kgds(1) = 0
       kgds(2) = in
       kgds(3) = jn
       kgds(5) = 0                             !first longitude 
       kgds(8) = int(360000./real(in)*(in-1))  !end longitude
       kgds(9) = int(360000./real(in))         !di longitudinal increment
       kgds(10) = int(180000./real(jn-1))      !dj latitudinal increment
       kgds(20) = 255    
       if (yrev) then
          kgds(4) = 90000                         !first latitude
          kgds(7) = -90000                        !end latitude
       else
          kgds(4) =-90000                        !first latitude
          kgds(7) = 90000                        !end latitude
          kgds(11) = 64                        !scanning mode flag
       end if
    else if (grid == 'gg') then !regular gaussian grid
       call gausslat(jn,slat,wlat)
       kgds(1) = 4
       kgds(2) = in
       kgds(3) = jn
       kgds(5) = 0                             !first longitude 
       kgds(8) = int(360000./real(in)*(in-1))  !end longitude
       kgds(9) = int(360000./real(in))         !di longitudinal increment
       kgds(10) = jn/2  !number of parallels between a pole and the equator
       kgds(20) = 255    
       if (yrev) then !N->S case
          kgds(4) = int(asin(slat(1))*180./pi*1000)
          kgds(7) = int(asin(slat(jn))*180./pi*1000)
       else           !S->N case
          kgds(4) = int(asin(slat(jn))*180./pi*1000)
          kgds(7) = int(asin(slat(1))*180./pi*1000)
          kgds(11) = 64                        !scanning mode flag
       end if
    else
       write(*,*) 'grid definition type is incorrect. grid=',grid
    end if

  end subroutine makekgds

end module ip_module
