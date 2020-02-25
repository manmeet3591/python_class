module read_module
  use type_module
  use param_module
  use constant_module, only: air_rd, earth_gravity
  use sigio_module, only: sigio_modpr
  use ip_module
  implicit none

  private

  public :: read_vcoord,read_grib1_mlev,read_grib1_plev, &
       &    read_grads,read_grib1_orography, read_grib2_plev, &
       &    read_grads_orography, read_o3clim

  interface getgrib1
     module procedure getgrib1_2d, getgrib1_3d
  end interface getgrib1
  interface getgrib2
     module procedure getgrib2_2d, getgrib2_3d
  end interface getgrib2

  integer,parameter :: nvarids = 14, max_char_size=200

contains

  subroutine read_vcoord(luvc,fname,vcoord)
    integer(kind=i4b),intent(in)  :: luvc
    character(len=*), intent(in)  :: fname
    real(kind=sp),    intent(out) :: vcoord(:,:)
    integer(kind=i4b)             :: n,k,iret
    
    iret = 0

    open(luvc,file=fname,iostat=iret)
    if (iret /= 0) then
       write(*,*) 'vertical level definition file open error! iret=',iret
       call abort
    end if
    read(luvc,'()') !skip header
    read(luvc,*,iostat=iret) ((vcoord(k,n),n=1,size(vcoord,2)),k=1,size(vcoord,1))
    if (iret /=0 ) then
       write(*,*) 'vertical level definition file read error! iret=',iret
       call abort
    end if
    close(luvc)
    
  end subroutine read_vcoord

  subroutine set_grib1_table(lutb, tablefile, varids)
    !!-----------------------------------------------------------------
    !! grib1のテーブルファイルを読み込んで変数のidを読み込む
    integer(kind=i4b), intent(in)  :: lutb
    character(len=*),  intent(in)  :: tablefile
    integer(kind=i4b), intent(out) :: varids(nvarids)
    
    integer(kind=i4b)            :: iret, varid
    character(len=5)             :: var
    character(len=max_char_size) :: line
     
    varids = 0
    
    open(lutb, file=tablefile, iostat=iret, status='old')
    if (iret/=0) then
       write(*,*) 'cannot find table definition file'
       call abort
    end if
    
    write(*,*)
    read(lutb,'(a)', iostat=iret) line
    write(*,*) trim(line)
    read(lutb,'(a)', iostat=iret) line
    write(*,*) trim(line)
    read(lutb,'(a)', iostat=iret) line
    write(*,*)
    
    varids = 0
    iret = 0
    do while(iret==0)
       read(lutb, *,iostat=iret) varid, var
       select case(trim(var))
       case ('zz') 
          varids(1) = varid
       case ('gh') 
          varids(2) = varid
       case ('tt') 
          varids(3) = varid
       case ('uu') 
          varids(4) = varid
       case ('vv') 
          varids(5) = varid
       case ('qq') 
          varids(6) = varid
       case ('rh') 
          varids(7) = varid
       case ('o3') 
          varids(8) = varid
       case ('cwc') 
          varids(9) = varid
       case ('clwc') 
          varids(10) = varid
       case ('ciwc') 
          varids(11) = varid
       case ('ps') 
          varids(12) = varid
       case ('lnsp') 
          varids(13) = varid
       case ('hs') 
          varids(14) = varid
       end select
    end do
    
    if (iret>0) then
       write(*,*) 'file read error. iret=', iret
       call abort
    end if
    
    close(lutb)
    print*, varids
  end subroutine set_grib1_table

  subroutine read_grib1_plev(lugb,fname,tablefile,yr,mn,dy,hr,in,jn,kn,levels,qvar,field,ps,kgdso)  
    integer(kind=i4b), intent(in)  :: lugb                ! file unit
    character(len=*),  intent(in)  :: fname, tablefile    ! file name
    integer(kind=i4b), intent(in)  :: yr,mn,dy,hr         ! year, month, day, hour
    integer(kind=i4b), intent(in)  :: in,jn,kn            ! x,y,z dimension size
    real(kind=i4b), intent(in)     :: levels(kn)          
    character(len=3),intent(in)    :: qvar
    real(kind=sp),intent(out)      :: field(in,jn,kn,7), ps(in,jn)
    integer(kind=i4b),intent(out)  :: kgdso(200)
    
    integer(kind=i4b) :: levs(kn), varids(nvarids)
    integer(kind=i4b) :: varid, levid
    real(kind=sp)     :: bufr(in,jn,kn)
    
    integer(kind=i4b) :: k, iret
    
    call set_grib1_table(lugb, tablefile, varids)

    call baopenr(lugb,fname,iret)
    if (iret /= 0) then
       write(*,*) 'file open error. iret=', iret
       call abort
    end if
    
    levid=100 !pressure level
    levs = int(levels)
    field = 0.
    bufr  = 0.

    write(*,*) 'read geopotential heigt'
    varid = varids(2)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,1),iret)
    if (iret/=0) then
       write(*,*) 'read geopotential'
       varid = varids(1)
       call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,1),iret)
       field(:,:,:,1) = field(:,:,:,1)/earth_gravity
    end if
    if (iret/=0) call abort

    write(*,*) 'read temperature'
    varid=varids(3)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,2),iret,kgdso=kgdso)
    if (iret/=0) call abort
    
    write(*,*) 'read u verocity'
    varid = varids(4)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,3),iret)
    if (iret/=0) call abort
    
    write(*,*) 'read v verocity'
    varid = varids(5)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,4),iret)
    if (iret /= 0) call abort
    
    if (trim(qvar) == 'rh') then
       varid = varids(7)
       write(*,*) 'read relative humidity'
    else if (trim(qvar) == 'qq') then
       varid = varids(6)
       write(*,*) 'read specific humidity'
    else
       write(*,*) 'input error'
       call abort
    end if
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,5),iret)
    if (iret < 0 .and. iret /= -kn) then
       write(*,*) 'request variables not found on some levels, set to zero'    
    else if (iret /= 0) then
       call abort 
    end if
    
    write(*,*) 'read ozon mixing ratio'
    varid = varids(8)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,6),iret)
    if (iret == -kn) then
       write(*,*) 'request variables not found, use climatological field'
    else if (iret < 0) then
       write(*,*) 'request variables not found on some levels, set to zero'
    else if (iret /= 0) then
       call abort 
    end if
    
    write(*,*) 'read cloud water mixing ratio'
    varid = varids(9)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,7),iret)
    if (iret < 0) then
       write(*,*) 'read cloud ice water mixing ratio'
       varid=varids(11)
       call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,7),iret)     
       write(*,*) 'read cloud liquid water mixing ratio'
       varid=varids(10)
       call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,bufr(:,:,:),iret)
       field(:,:,:,7) = field(:,:,:,7) + bufr(:,:,:)
    end if
    if (iret < 0) then
       write(*,*) 'request variables not found on some levels, set to zero'    
    else if (iret /= 0) then
       call abort 
    end if

    write(*,*) 'read surface pressure'
    varid=varids(12)
    call getgrib1(lugb,in,jn,-1,varid,-1,yr,mn,dy,hr,ps(:,:),iret)
    if (iret /= 0) then
       write(*,*) 'read log surface pressure'
       varid=varids(13)
       call getgrib1(lugb,in,jn,-1,varid,-1,yr,mn,dy,hr,ps(:,:),iret)
       ps = exp(ps)
    end if
    if (iret/=0) then
       write(*,*) 'surface presure not found'
       ps = -9999.
    end if
    
    call baclose(lugb,iret)

  end subroutine read_grib1_plev

  subroutine read_grib1_mlev(lugb,fname,tablefile,yr,mn,dy,hr,in,jn,kn,levels,qvar,  &
       &                     zvar,field,ps,kgdso)  
    !!--------------------------------------------------------------------------------------
    !! grib1ファイルを読み込む（モデル面）
    !!--------------------------------------------------------------------------------------
    integer(kind=i4b), intent(in)  :: lugb
    character(len=*),  intent(in)  :: fname, tablefile
    integer(kind=i4b), intent(in)  :: yr,mn,dy,hr
    integer(kind=i4b), intent(in)  :: in,jn,kn
    real(kind=i4b), intent(in)     :: levels(kn)
    character(len=2),intent(out)   :: zvar
    character(len=3),intent(in)    :: qvar
    real(kind=sp),intent(out)      :: field(in,jn,kn,7),ps(in,jn)
    integer(kind=i4b),intent(out)  :: kgdso(200)
    
    real(kind=sp)     :: gg = air_rd/earth_gravity
    integer(kind=i4b) :: levs(kn),varids(nvarids)
    integer(kind=i4b) :: varid, levid
    real(kind=sp)     :: bufr(in,jn,kn)
    
    integer(kind=i4b) :: k, iret
    
    call set_grib1_table(lugb,tablefile,varids)
    call baopenr(lugb,fname,iret)
    if (iret /= 0) then
       write(*,*) 'file open error. iret=',iret
       call abort
    end if
    
    levid=109 !hybrid level
    levs = int(levels)
    field = 0.
    ps = 0.
    bufr  = 0.
    
    write(*,*) 'read surface pressure'
    varid=varids(12)
    call getgrib1(lugb,in,jn,-1,varid,-1,yr,mn,dy,hr,ps(:,:),iret)
    if (iret /= 0) then
       write(*,*) 'read log surface pressure'
       varid=varids(13)
       call getgrib1(lugb,in,jn,-1,varid,-1,yr,mn,dy,hr,ps(:,:),iret)
       ps = exp(ps)
    end if
    if (iret/=0) call abort
    
    write(*,*) 'read geopotential heigt'
    varid = varids(2)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,1),iret)
    zvar = 'gh'
    if (iret /= 0) then
       write(*,*) 'read geopotential'
       varid = varids(1)
       call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,1),iret)
       field(:,:,:,1) = field(:,:,:,1)/earth_gravity
       zvar = 'gh'
    end if
    if (iret /= 0) then
       write(*,*) 'read lowest geopotential'
       varid = varids(1)
       call getgrib1(lugb,in,jn,1,varid,levid,yr,mn,dy,hr,field(:,:,1,1),iret)
       field(:,:,1,1) = field(:,:,1,1)/earth_gravity
       zvar = 'z1' 
    end if
    if (iret /= 0 ) then
       write(*,*) 'request field not found. caluclate geopotential height using surface height'
       write(*,*) 'read geometric height'
       varid = varids(14)
       call getgrib1(lugb,in,jn,1,varid,1,yr,mn,dy,hr,field(:,:,1,1),iret)
       field(:,:,1,1) = field(:,:,1,1)/earth_gravity
       zvar = 'hs'
    end if

    write(*,*) 'read temperature'
    varid=varids(3)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,2),iret,kgdso=kgdso)
    if (iret/=0) call abort
    
    write(*,*) 'read u verocity'
    varid=varids(4)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,3),iret)
    if (iret/=0) call abort
    
    write(*,*) 'read v verocity'
    varid=varids(5)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,4),iret)
    if (iret/=0) call abort
    
    if (trim(qvar) == 'rh') then
       varid = varids(7)
       write(*,*) 'read relative humidity'
    else if (trim(qvar) == 'qq') then
       varid = varids(6)
       write(*,*) 'read specific humidity'
    else
       write(*,*) 'input error'
       call abort
    end if
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,5),iret)
    if (iret < 0 .and. iret /= -kn) then
       write(*,*) 'request variables not found on some levels, set to zero'    
    else if (iret /= 0) then
       call abort 
    end if
    
    write(*,*) 'read ozon mixing ratio'
    varid = varids(8)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,6),iret)
    if (iret == -kn) then
       write(*,*) 'request variables not found, use climatological field'
    else if (iret < 0) then
       write(*,*) 'request variables not found on some levels, set to zero'
    else if (iret /= 0) then
       call abort 
    end if
    
    write(*,*) 'read cloud water mixing ratio'
    varid = varids(9)
    call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,7),iret)
    if (iret == -kn) then
       write(*,*) 'read cloud ice water mixing ratio'
       varid=varids(11)
       call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,field(:,:,:,7),iret)     
       write(*,*) 'read cloud liquid water mixing ratio'
       varid=varids(10)
       call getgrib1(lugb,in,jn,kn,levs,varid,levid,yr,mn,dy,hr,bufr(:,:,:),iret)
       field(:,:,:,7) = field(:,:,:,7) + bufr(:,:,:)
    end if
    if (iret < 0) then
       write(*,*) 'request variables not found on some levels, set to zero'    
    else if (iret /= 0) then
       call abort 
    end if
    
    call baclose(lugb,iret)

  end subroutine read_grib1_mlev

  subroutine read_grib2_plev(lugb,fname,yr,mn,dy,hr,in,jn,kn,levels,qvar,field,ps,kgdso)
    integer(kind=i4b), intent(in) :: lugb
    character(len=*), intent(in)  :: fname
    integer(kind=i4b), intent(in) :: yr,mn,dy,hr
    integer(kind=i4b), intent(in) :: in,jn,kn
    real(kind=i4b), intent(in)    :: levels(kn)
    real(kind=sp),intent(out)     :: field(in,jn,kn,7), ps(in,jn)
    character(len=2),intent(in)   :: qvar
    integer(kind=i4b),intent(out) :: kgdso(200)
    integer(kind=i4b) :: levs(kn)
    integer(kind=i4b) :: varid1,varid2,levid
    integer(kind=i4b) :: iret
    logical :: pass99 = .true. ! test
    
    !pressure level (hPa)
    levs = int(levels)
    
    call baopenr(lugb,fname,iret)
    if (iret /= 0) then
       write(*,*) 'file open error. iret=',iret
       call abort
    end if
    !initialize
    field = 0.
    
    !read field on layer
    levid=100 !pressure level
    write(*,*) 'read geopotential height'
    varid1=3
    varid2=5
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,1),iret,kgdso=kgdso)  
    if (iret/=0) call abort

    write(*,*) 'read temperature'
    varid1=0
    varid2=0
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,2),iret)
    if (iret/=0) call abort

    write(*,*) 'read u verocity'
    varid1=2
    varid2=2
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,3),iret)
    if (iret/=0) call abort

    write(*,*) 'read v verocity'
    varid1=2
    varid2=3
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,4),iret)
    if (iret/=0) call abort

    if (qvar=='rh') then
       write(*,*) 'read relative humidity'
       varid1=1
       varid2=1
    else if (qvar=='qq') then
       write(*,*) 'read specific humidity'
       varid1=1
       varid2=0
    else
       write(*,*) 'input error'
       call abort
    end if
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,5),iret)
    if (iret == 99) then
       write(*,*) 'request variables not found on some levels, set to zero'
    else if (iret /= 0) then
       call abort 
    end if

    write(*,*) 'read ozone mixing ratio'
    varid1=14
    varid2=1 
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,6),iret)
    if (iret == 99) then
       write(*,*) 'request variables not found on some levels, set to zero'
    else if (iret == 999) then
       write(*,*) 'request variables not found, use climatological field'
    else if (iret /= 0) then
       call abort 
    end if

    write(*,*) 'read cloud water mixing ratio'
    varid1=1
    varid2=22
    call getgrib2(lugb,in,jn,kn,levs,varid1,varid2,levid,yr,mn,dy,hr,field(:,:,:,7),iret)
    if (iret == 99) then
       write(*,*) 'request variables not found on some levels, set to zero'
    else if (iret /= 0) then
       call abort 
    end if

    write(*,*) 'read surface pressure'
    varid1=3
    varid2=0
    levid =1  !surface level
    call getgrib2(lugb,in,jn,-1,varid1,varid2,levid,yr,mn,dy,hr,ps,iret)
    if (iret /= 0) then
       write(*,*) 'surface pressure not found'
       ps = -9999.
    end if

    call baclose(lugb,iret)
  
end subroutine read_grib2_plev

  subroutine read_grads(lugb,fname,in,jn,kn,recs,psrec,field,ps,big_endian)
    integer(kind=i4b), intent(in) :: lugb
    character(len=*), intent(in)  :: fname
    integer(kind=i4b), intent(in) :: in,jn,kn
    integer(kind=i4b), intent(in) :: recs(7), psrec
    real(kind=sp),intent(out)     :: field(in,jn,kn,7), ps(in,jn)
    logical, optional             :: big_endian
    integer(kind=i4b) :: iret,i,rec,k

    if (present(big_endian) .and. big_endian) then
       open(lugb,file=fname,access='direct',form='unformatted',status='old', &
            & iostat=iret,recl=4*in*jn,convert='big_endian')
    else
       open(lugb,file=fname,access='direct',form='unformatted',status='old', &
            & iostat=iret,recl=4*in*jn)
    end if
    if (iret /= 0) then
       write(*,*) 'file open error. iret=',iret
       call abort
    end if
    field = 0.
    do i = 1, 7
       if (recs(i) < 0) then
          write(*,*) 'no field. set to zero.'
       else
          rec = recs(i)
          if (rec == 0) rec = (i-1)*kn+1
          write(*,*) 'read rec=',rec
          do k = 1, kn
             read(lugb,rec=rec+k-1) field(:,:,k,i)
          end do
       end if
    end do

    ps = 0.
    if (psrec >= 0) then
       rec=psrec
       if (rec == 0) rec=7*kn+1       
       read(lugb,rec=rec) ps(:,:)
    end if
    close(lugb)

  end subroutine read_grads

  subroutine read_grib1_orography(lugb,fname,in,jn,data,kgdso)
    integer(kind=i4b), intent(in) :: lugb
    character(len=*),  intent(in) :: fname
    integer(kind=i4b), intent(in) :: in,jn
    real(kind=sp),intent(out)     :: data(in,jn)
    integer(kind=i4b),intent(out),optional :: kgdso(200)
    integer(kind=i4b) :: iret,kgds(200),varids(nvarids),varid,levid,grib_version
    
    call baopenr(lugb,fname,iret)
    if (iret /= 0) then
       write(*,*) 'orography file open error'
       call abort
    end if

    levid=1
    write(*,*) 'read surface geometric height'

    varid=8 !geometric height (NCEP)
    kgds=0
    data=0.
    call getgrib1(lugb,in,jn,-1,varid,levid,-1,-1,-1,-1,data,iret,kgdso=kgds)
    if (iret /= 0) then
       write(*,*) 'read error. iret=', iret
       call abort
    end if
    call baclose(lugb,iret)
    if (present(kgdso)) then
       kgdso = kgds
    end if

  end subroutine read_grib1_orography

  subroutine read_grads_orography(lugb,fname,in,jn,hrec,field,big_endian)
    integer(kind=i4b), intent(in) :: lugb
    character(len=*), intent(in)  :: fname
    integer(kind=i4b), intent(in) :: in,jn,hrec
    real(kind=sp),intent(out)     :: field(in,jn)
    logical, optional             :: big_endian
    integer(kind=i4b) :: iret,rec
    
    if (present(big_endian) .and. big_endian) then
       open(lugb,file=fname,access='direct',form='unformatted',status='old', &
            & iostat=iret,recl=4*in*jn,convert='big_endian')
    else
       open(lugb,file=fname,access='direct',form='unformatted',status='old', &
            & iostat=iret,recl=4*in*jn)
    end if
    if (iret /= 0) then
       write(*,*) 'file open error. iret=',iret
    end if
    field = 0.
    if (rec==0) then
       rec=1
    else
       rec=hrec
    end if
    write(*,*) 'read rec=',rec
    read(lugb,rec=rec) field(:,:)
    close(lugb)
    
  end subroutine read_grads_orography

  subroutine getgrib1_3d(lugb,in,jn,kn,levels,varid,levid,yr,mn,dy,hr,odata,iret,kgdso)
    integer(kind=i4b), intent(in) :: lugb
    integer(kind=i4b), intent(in) :: in, jn, kn
    integer(kind=i4b), intent(in) :: levels(kn)
    integer(kind=i4b), intent(in) :: varid,levid,yr,mn,dy,hr
    integer(kind=i4b), intent(out):: iret
    real(kind=sp),intent(out)     :: odata(in,jn,kn)
    integer,intent(out),optional  :: kgdso(200)
    
    integer(kind=i4b) i,j,k,nfiret
    integer(kind=i4b) lugi,jskip,jpds(200),jgds(200),ndata,lskip,kpds(200),kgds(200)
    logical lb(in*jn)
    real(kind=sp) temp(jn)
    
    lugi=0  !no index file
    jskip=0 !serach for beginning
    jpds=-1
    jgds=-1
    jpds(5)=varid
    jpds(6)=levid
    if (yr == -1) then
       jpds(8)=-1
    else
       jpds(8)=mod(yr,100) !year of century
    end if
    jpds(9)=mn
    jpds(10)=dy
    jpds(11)=hr
    iret = 0
    nfiret = 0             ! request not found at all level flat
    do k = 1, kn
       jpds(7)=levels(k)
       ndata=0
       kpds=0
       kgds=0     
       call getgb(lugb,lugi,in*jn,jskip,jpds,jgds,ndata,lskip,kpds,kgds,lb,odata(:,:,k),iret)
       if (iret == 99) then
          write(*,*) 'request not found.'
          write(*,*) 'jpds(5)=',jpds(5),',jpds(6)=',jpds(6),',jpds(7)=',jpds(7)
          nfiret = nfiret + 1
          cycle
       else if (iret /= 0) then
          write(*,*) 'file read error. iret=',iret          
          return
       end if
       if (ndata == 0) then
          write(*,*) 'error: no data field on levs=',jpds(7) 
          iret = 101
          return
       end if
       if (kgds(2) /= in .or. kgds(3) /= jn) then
          write(*,*) 'file field size is incorrect.'
          write(*,*) 'request field size is xn=', in, 'yn=', jn
          write(*,*) 'file field size is xn=', kgds(2), 'yn=', kgds(3)
          iret = 102
          return
       end if
    end do
    if (present(kgdso)) then
       kgdso=kgds
    end if

    ! number of not found request is returned as iret 
    iret = -nfiret     
    
  end subroutine getgrib1_3d
  
  subroutine getgrib1_2d(lugb,in,jn,level,varid,levid,yr,mn,dy,hr,odata,iret,kgdso)
    integer(kind=i4b), intent(in) :: lugb
    integer(kind=i4b), intent(in) :: in, jn, level
    integer(kind=i4b), intent(in) :: varid,levid,yr,mn,dy,hr
    real(kind=sp),    intent(out) :: odata(in,jn)
    integer(kind=i4b),intent(out),optional :: kgdso(200)
    integer(kind=i4b), intent(out):: iret  
    integer(kind=i4b) :: i,j
    integer(kind=i4b) :: lugi,jskip,jpds(200),jgds(200),ndata,lskip,kpds(200),kgds(200)
    logical :: lb(in*jn)
    real(kind=sp) :: temp(jn)
    
    lugi=0  !no index file
    jskip=0 !serach for beginning
    jpds=-1
    jgds=-1
    jpds(5)=varid
    jpds(6)=levid
    jpds(7)=level      !(-1=wild card)
    if (yr == -1) then
       jpds(8)=-1
    else
       jpds(8)=mod(yr,100) !year of century
    end if
    jpds(9)=mn
    jpds(10)=dy
    jpds(11)=hr
    iret = 0
    kpds = 0
    kgds = 0
    odata = 0.
    call getgb(lugb,lugi,in*jn,jskip,jpds,jgds,ndata,lskip,kpds,kgds,lb,odata,iret)
    if (iret /= 0) then
       if (iret == 99) then
          write(*,*) 'request not found. jpds(5)=',jpds(5),',jpds(6)=',jpds(6),',jpds(7)=',jpds(7)
          return
       else
          write(*,*) 'file read error. iret=',iret
       end if
    end if
    if (ndata == 0) then
       write(*,*) 'error: no data field.'
       iret = 101
    end if
    if (kgds(2) /= in .or. kgds(3) /= jn) then
       write(*,*) 'file field size is incorrect.'
       write(*,*) 'request field size is xn=', in, 'yn=', jn
       write(*,*) 'file field size is xn=', kgds(2), 'yn=', kgds(3)
       iret = 102
    end if
    if (present(kgdso)) then
       kgdso=kgds
    end if
  end subroutine getgrib1_2d
  
  subroutine getgrib2_3d(lugb,in,jn,kn,levels,varid1,varid2,levid,yr,mn,dy,hr,odata,iret,kgdso)
    use grib_mod
    integer(kind=i4b), intent(in) :: lugb
    integer(kind=i4b), intent(in) :: in, jn, kn
    integer(kind=i4b), intent(in) :: levels(kn)
    integer(kind=i4b), intent(in) :: varid1,varid2,levid,yr,mn,dy,hr
    real(kind=sp),intent(out)     :: odata(in,jn,kn)
    integer(kind=i4b),intent(out),optional :: kgdso(200)
    
    integer(kind=i4b) i,j,k,iret,n,nfiret,igds(5),igrid
    integer(kind=i4b) lugi,jskip,jdisc,jids(200),jpdtn,jpdt(200),jgdtn,jgdt(200),ndata
    logical unpack
    type(gribfield) gfld
    real(kind=sp) temp(jn)
    
    lugi=0     !no index file
    jskip=0    !serach for beginning
    jdisc=-1   !description number
    jids=-9999
    jpdtn=0
    jpdt=-9999
    jgdtn=-1
    jgdt=-9999
    unpack=.true.
    
    jids(6)=yr
    jids(7)=mn
    jids(8)=dy
    jids(9)=hr
    
    jpdt(1)=varid1
    jpdt(2)=varid2
    jpdt(10)=levid
    !  jpdt(11)=0
    iret = 0
    nfiret = 0
    do k = 1, kn
       jpdt(12)=levels(k)
       call getgb2(lugb,lugi,jskip,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,unpack,ndata,gfld,iret)       
       if (iret == 99) then
          write(*,*) 'request not found.  var=',jpdt(1),jpdt(2),',level=',jpdt(12)
          write(*,*) 'no field. set to zero.'
          cycle   
       else if (iret /= 0) then
          write(*,*) 'file read error. iret=',iret
          return
       end if
       if (ndata == 0) then
          write(*,*) 'error: no data field on levs=',jpdt(12)        
          iret = 101
          return
       end if
       if (gfld%ngrdpts /= in*jn) then
          write(*,*) 'file field size is incorrect.'
          write(*,*) 'request field size is xn*yn=',in*jn
          write(*,*) 'file field size is xn*yn=', gfld%ngrdpts
          iret = 102
          return
       end if
       n=1
       do i = 1, in
          do j = 1, jn
             odata(i,j,k) = gfld%fld(n)
             n=n+1
          end do
       end do
       if (gfld%igdtnum/=0 .and. gfld%igdtnum/=4) then
          write(*,*) 'This file grid neigher lon/lat or gaussian. igdtnum=',gfld%igdtnum
          iret = 103
          return
       end if
       if (present(kgdso)) then
          !see more detail in gdt2gds.f and getgb2.f of g2lib
          igds(1) = gfld%griddef
          igds(2) = gfld%ngrdpts
          igds(3) = gfld%numoct_opt
          igds(4) = gfld%interp_opt
          igds(5) = gfld%igdtnum
          call gdt2gds(igds,gfld%igdtmpl,gfld%num_opt,gfld%list_opt,kgdso,igrid,iret)
       end if
    end do

    ! number of not found request is returned as iret 
    iret = -nfiret     
    
  end subroutine getgrib2_3d
  subroutine getgrib2_2d(lugb,in,jn,level,varid1,varid2,levid,yr,mn,dy,hr,odata,iret,kgdso)
    use grib_mod
    integer(kind=i4b), intent(in) :: lugb
    integer(kind=i4b), intent(in) :: in, jn
    integer(kind=i4b), intent(in) :: level
    integer(kind=i4b), intent(in) :: varid1,varid2,levid,yr,mn,dy,hr
    real(kind=sp),intent(out)     :: odata(in,jn)
    integer(kind=i4b),intent(out),optional :: kgdso(200)
    
    integer(kind=i4b) i,j,iret,n,igds(5),igrid
    integer(kind=i4b) lugi,jskip,jdisc,jids(200),jpdtn,jpdt(200),jgdtn,jgdt(200),ndata
    logical unpack
    type(gribfield) gfld
    real(kind=sp) temp(jn)
    
    lugi=0     !no index file
    jskip=0    !serach for beginning
    jdisc=-1   !description number
    jids=-9999
    jpdtn=0
    jpdt=-9999
    jgdtn=-1
    jgdt=-9999
    unpack=.true.
    
    jids(6)=yr
    jids(7)=mn
    jids(8)=dy
    jids(9)=hr
    
    jpdt(1)=varid1
    jpdt(2)=varid2
    jpdt(10)=levid
    !  jpdt(11)=0
    jpdt(12)=level
    iret = 0

    call getgb2(lugb,lugi,jskip,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,unpack,ndata,gfld,iret)
    if (iret == 99) then
       write(*,*) 'request not found.  var=',jpdt(1),jpdt(2),',level=',jpdt(12)
       return
    else if (iret /= 0) then
       write(*,*) 'file read error. iret=',iret
       return
    end if
    if (ndata == 0) then
       write(*,*) 'error: no data field on levs=',jpdt(12)        
       iret = 101
       return
    end if
    if (gfld%ngrdpts /= in*jn) then
       write(*,*) 'file field size is incorrect.'
       write(*,*) 'request field size is xn*yn=',in*jn
       write(*,*) 'file field size is xn*yn=', gfld%ngrdpts
       iret = 102
       return
    end if
    n=1
    do i = 1, in
       do j = 1, jn
          odata(i,j) = gfld%fld(n)
          n=n+1
       end do
    end do
    if (gfld%igdtnum/=0 .and. gfld%igdtnum/=4) then
       write(*,*) 'This file grid neigher lon/lat or gaussian. igdtnum=',gfld%igdtnum
       iret = 103
       return
    end if
    if (present(kgdso)) then
       !see more detail in gdt2gds.f and getgb2.f of g2lib
       igds(1) = gfld%griddef
       igds(2) = gfld%ngrdpts
       igds(3) = gfld%numoct_opt
       igds(4) = gfld%interp_opt
       igds(5) = gfld%igdtnum
       call gdt2gds(igds,gfld%igdtmpl,gfld%num_opt,gfld%list_opt,kgdso,igrid,iret)
    end if

  end subroutine getgrib2_2d


  subroutine read_o3clim(luo3,fname,yr,mn,dy,hr,kgdso,o3out,plevs)
    !copy and modified from GETO3C in chgres.f
    integer(kind=i4b), parameter  :: kmc=17, jmc=18 !ozone clim fix file levels, lats
    integer(kind=i4b), intent(in) :: luo3,yr,mn,dy,hr
    character(len=*),intent(in)   :: fname
    integer(kind=i4b),intent(in)  :: kgdso(200)
    real(kind=sp),intent(out)     :: o3out(kgdso(2),kgdso(3),kmc), plevs(kgdso(2),kgdso(3),kmc)
    
    integer(kind=i4b) :: i,j,k,mon,iret,mn1,mn2,dy1,dy2
    integer(kind=i4b) :: mndys(12), kgdsi(200)
    real(kind=sp)     :: lats(jmc),bufr(jmc,kmc,12),o3clim(2,jmc,kmc),plevc(kmc)
    
    data mndys /31,28,31,30,31,30,31,31,30,31,30,31/
    
    open(luo3,file=fname,status='old',iostat=iret)
    if (iret /= 0) then
       write(*,*) 'file open error. iret=',iret
    end if
    ! read clim pressure
    read(luo3,*) plevc
    do i = 1, 12
       read(luo3,*) (mon,lats(j),(bufr(j,k,mon),k=1,kmc),j=1,jmc)
    end do
    close(luo3)
    !convert from PPM to Kg/Kg
    bufr = bufr*1.655e-6

    !reverse
    do k = 1, kmc
       plevs(:,:,k) = plevc(kmc-k+1)*100 !hPa->Pa
       bufr(:,k,:) = bufr(:,kmc-k+1,:)
    end do

    !time interpolate
    !date in first half of month; interpolate with the previous month
    if (dy<=mndys(mn)/2) then
       if (mn==1) then
          mn1 = 12
          mn2 = 1
       else
          mn1 = mn-1
          mn2 = mn
       end if
       dy2=mndys(mn2)/2-dy
       dy1=mndys(mn1)-mndys(mn1)/2+dy
    else
       if (mn==12) then
          mn1=12
          mn2=1
       else
          mn1=mn
          mn2=mn+1
       end if
       dy1=dy-mndys(mn1)/2.
       dy2=mndys(mn2)/2+mndys(mn1)-dy
    end if
    o3clim(1,:,:)=(dy2*bufr(:,:,mn1)+dy1*bufr(:,:,mn2))/(dy1+dy2)
    o3clim(2,:,:)=(dy2*bufr(:,:,mn1)+dy1*bufr(:,:,mn2))/(dy1+dy2)

    !interpolate
    kgdsi=0
    kgdsi(2)=2
    kgdsi(3)=jmc
    kgdsi(4)=int(lats(1)*1000)
    kgdsi(5)=0
    kgdsi(6)=int(lats(jmc)*1000)
    kgdsi(7)=180000
    kgdsi(11)=64
    kgdsi(20)=255
    call kgds2kgds(0,o3clim,kgdsi,kgdso,o3out)

  end subroutine read_o3clim

  
end module read_module
