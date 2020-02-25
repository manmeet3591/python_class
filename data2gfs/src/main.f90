program main
  use type_module
  use param_module
  use read_module
  use vcoord_def_module, only: vcoord=>gfs_hyblev_l64
  use sigio_module
  use sigio_r_module
  use calcmet_module
  use ip_module
  use p2p_module
  use gfs_module
  use constant_module, only: air_rd, earth_gravity
  implicit none

  integer(kind=i4b)              :: in, jn, kn, yr, mn, dy, hr, recs(7), idusr,psrec, hrec
  real(kind=sp)                  :: fhour
  real(kind=sp), allocatable     :: levels(:)
  character(len=5)               :: form
  character(len=2)               :: grid, levtyp
  logical           :: o3clim, yrev, big_endian
  real(kind=sp)     :: undef

  integer(kind=i4b)              :: jcap, lonr, latr, levs,maxwv,iptype
  
  integer(kind=i4b)            :: kgdso(200),kgdsi(200)
  real(kind=sp),allocatable    :: bufr(:,:,:,:),plbufr(:,:,:),fldin(:,:,:,:),plin(:,:,:),vcoordin(:,:)
  character(len=3)             :: qvar
  character(len=2)             :: zvar
  real(kind=sp),allocatable    :: fldout(:,:,:,:),plout(:,:,:),hs(:,:),ps(:,:)
  type(sigio_head)  :: head
  type(sigio_data)  :: data

  integer(kind=i4b) :: i,k, iret

  namelist/gfsparam/jcap,lonr,latr,levs,idusr,o3clim
  namelist/dataparam/in,jn,kn,yr,mn,dy,hr,fhour,levtyp,form,qvar,yrev,recs,big_endian,undef,grid,psrec,hrec,zvar
  namelist/ipparam/iptype,maxwv
  namelist/level/levels

  !--------------------------------------
  ! read namelist
  !
  qvar='qq '
  zvar='gh'
  levtyp='pl'
  o3clim=.false.
  fhour=0.
  ! for interpolation parameter
  iptype=4
  maxwv=-1
  ! for binary mode parameter
  grid='ll'           !lat-lon grid
  yrev=.false.        !south to north
  big_endian=.false.
  recs=0
  psrec=0
  hrec=0
  undef=9.999e+20

  write(*,*) 'read namelist'
  open(lunm,file='data2gfs_namelist')
  read(lunm,nml=gfsparam)
  read(lunm,nml=dataparam)
  read(lunm,nml=ipparam)
  allocate(levels(kn))
  read(lunm,nml=level)
  close(lunm)
  write(*,*) 'input data'
  write(*,'(A8,A5,A5,A5)') ' levtyp=',levtyp,' form=',form
  write(*,'(A4,I4,A4,I4,A4,I4)') ' xn=',in,' yn=',jn,' zn=',kn
  write(*,'(A6,I4,A7,I4,A5,I4,A6,I4)') ' year=',yr,' month=',mn,' day=',dy,' hour=',hr
  write(*,*)
  write(*,*) 'output'
  write(*,'(A6,I4,A6,I4,A6,I4,A6,I3)') ' jcap=',jcap,' lonr=',lonr,' latr=',latr,' levs=',levs
  write(*,*)
  write(*,*) 'PASS: read namelist'
  write(*,*) 

  
  !------------------------------
  ! read input field

  write(*,*) 'read input file ...'
  allocate(bufr(in,jn,kn,7),ps(in,jn))
  allocate(fldin(lonr,latr,kn,7),plin(lonr,latr,kn))
  bufr  = 0.
  ps    = 0.
  plin  = 0.
  fldin = 0.

  if (form == 'grib1') then
     if (levtyp == 'hl') then
        call read_grib1_mlev(lugb,'input.grib','input.grib1.table',yr,mn,dy,hr,in,jn,kn,levels,qvar, &
             &               zvar,bufr,ps,kgdso=kgdsi)
     else if (levtyp == 'pl') then
        call read_grib1_plev(lugb,'input.grib','input.grib1.table',yr,mn,dy,hr,in,jn,kn,levels,qvar, &
             &               bufr,ps,kgdso=kgdsi)
     else
        write(*,*) 'levtyp keyword is incorrect'
        call abort
     end if
  else if (form == 'grib2') then
     if (levtyp == 'pl') then
        call read_grib2_plev(lugb,'input.grib2',yr,mn,dy,hr,in,jn,kn,levels,qvar,bufr,ps,kgdso=kgdsi)
     else
        write(*,*) 'levtyp keyword is incorrect'
        call abort
     end if
  else if (form == 'grads') then
     if (levtyp == 'hl' .or. levtyp == 'sl' .or. levtyp == 'pl') then
        call read_grads(lugb,'input.bin',in,jn,kn,recs,psrec,bufr,ps,big_endian=big_endian)
     else 
        write(*,*) 'levtyp keyword is incorrect'
        call abort
     end if
     call makekgds(grid,in,jn,yrev,kgdsi)
  else
     write(*,*) 'this format is unabaialbe. form=',form
     call abort
  end if
  write(*,*) ' kgds(1:10)='
  write(*,*) (kgdsi(i),i=1,10)
  write(*,*) ' kgds(11:20)='
  write(*,*) (kgdsi(i),i=11,20)
  write(*,*)
  write(*,*) 'PASS: read input file'
  write(*,*) 
  where(abs(bufr)>=undef) bufr=0.  ! undef value force to set zero

  !-------------------------------------------------------------
  ! calc input file pressure level
  !
  write(*,*) 'calculate input file pressure level ...'
  if (levtyp == 'pl') then
     do k = 1, kn
        plin(:,:,k) = levels(k)*100.
     end do
  else if (levtyp == 'hl') then
     allocate(vcoordin(kn+1,2),plbufr(in,jn,kn))
     call read_vcoord(luvc,'input.vcoord.txt',vcoordin)
     call sigio_modpr(in*jn,in*jn,kn,2,2,2,vcoordin,iret, &
          &           ps=ps,pm=plbufr)
  else if (levtyp == 'sl') then
     allocate(vcoordin(kn,1),plbufr(in,jn,kn))
     call read_vcoord(luvc,'input.vcoord.txt',vcoordin)
     do k = 1, kn
        plbufr(:,:,k) = ps*vcoordin(k,1)
     end do
  end if
  write(*,*) 'PASS: calculate input file pressure level'
  write(*,*)
  
  !---------------------------------------------------
  ! read input file orography and calculate geopotential height if no geopotential case
  !
  if (zvar=='z1') then !input lowest level geopotential case
     do k = 2, kn
        bufr(:,:,k,1) = bufr(:,:,k-1,1) - 0.5*air_rd/earth_gravity*(bufr(:,:,k-1,2) &
             &           +bufr(:,:,k,2))*(log(plbufr(:,:,k))-log(plbufr(:,:,k-1)))
     end do
  else if (zvar=='hs') then !no input geopotential case
     allocate(hs(in,jn))
     write(*,*) 'calculate geopotential height from surface height. ...'
     if (form=='grib1') then
        hs(:,:) = bufr(:,:,1,1)
     else if (form=='grads') then
        call read_grads_orography(lugboi,'input.geo.bin',in,jn,hrec,hs,big_endian=big_endian)
     end if
     bufr(:,:,:,1) = calcmet_z(bufr(:,:,:,2),hs,ps,plbufr)     
     write(*,*) 'PASS: caculate geopotential height'
     write(*,*)
     deallocate(hs)
  end if

  !-----------------------------
  ! new orography
  !
  allocate(hs(lonr,latr))
  write(*,*) 'read new orography ...'
  call read_grib1_orography(61,'input.orography.grib',lonr,latr,hs,kgdso=kgdso)
  write(*,*) ' kgds(1:10)='
  write(*,*) (kgdso(i),i=1,10)
  write(*,*) ' kgds(11:20)='
  write(*,*) (kgdso(i),i=11,20)
  write(*,*)
  write(*,*) 'PASS: read new orography'
  write(*,*) 

  !-----------------------------
  ! horizontaly interpolate
  ! iptype = 0: bilinear 1:cubic 4:spectral
  write(*,*) 'horizontaly interpolate ...'
  call kgds2kgds(iptype,bufr(:,:,:,1),kgdsi,kgdso,fldin(:,:,:,1),maxwv=maxwv) !geopotential height
  call kgds2kgds(iptype,bufr(:,:,:,2),kgdsi,kgdso,fldin(:,:,:,2),maxwv=maxwv) !temperature
  call kgds2kgds(iptype,bufr(:,:,:,3),bufr(:,:,:,4),kgdsi,kgdso, &
       &         fldin(:,:,:,3),fldin(:,:,:,4),maxwv=maxwv)                   !u,v verocity
  call kgds2kgds(iptype,bufr(:,:,:,5),kgdsi,kgdso,fldin(:,:,:,5),maxwv=maxwv) !water vapor variable
  call kgds2kgds(iptype,bufr(:,:,:,6),kgdsi,kgdso,fldin(:,:,:,6),maxwv=maxwv) !ozone
  call kgds2kgds(iptype,bufr(:,:,:,7),kgdsi,kgdso,fldin(:,:,:,7),maxwv=maxwv) !cloud water ratio
  if (levtyp=='hl' .or. levtyp=='sl') then
     call kgds2kgds(iptype,plbufr,kgdsi,kgdso,plin(:,:,:),maxwv=maxwv)        !pressure
     deallocate(plbufr)
  end if
  write(*,*) 'PASS: horizontaly interpolate'
  write(*,*)
  deallocate(bufr,ps)


  !-----------------------------
  !calculate new surface pressure
  !
  allocate(ps(lonr,latr))
  write(*,*) 'calc new surface pressure ...'
  ps(:,:) = calcmet_ps(plin,fldin(:,:,:,2),fldin(:,:,:,1),hs)
  write(*,*) 'PASS: calc new surface pressure'
  write(*,*) 

  !----------------------------
  !calculate new pressure level
  !
  allocate(plout(lonr,latr,levs))
  write(*,*) 'calculate new model pressure ...'
  call sigio_modpr(lonr*latr,lonr*latr,levs,nvcoord,idvc,idsl,vcoord,iret, &
       &          ps=ps(:,:),pm=plout(:,:,:))
  write(*,*) 'PASS: calclate new model pressure'
  write(*,*) 

  !---------------------------------------------------------
  !read and interpolate climatorogical ozone field if need.
  if (o3clim) then
     allocate(bufr(lonr,latr,17,2)) !bufr(:,:,:,1):ozone field, bufr(:,:,:,2):pressure
     write(*,*) 'climatorogical ozone field input...'
     call read_o3clim(luo3,'input.o3clim.txt',yr,mn,dy,hr,kgdso,bufr(:,:,:,1),bufr(:,:,:,2))
     write(*,*) 'PASS: climatorogical ozone input'
     write(*,*)
  end if

  !---------------------------------------------------------------------
  !verticaly interpolate input pressure level -> output pressure level
  !
  write(*,*) 'vertical interpolation ...'
  allocate(fldout(lonr,latr,levs,6))
  fldout = 0.
  call p2p_extrapolate_T(fldin(:,:,:,2),plin,plout,ps,ps,hs,fldout(:,:,:,1))  !temperature
  call p2p(fldin(:,:,:,3),plin,plout,ps,fldout(:,:,:,2))           !u-wind
  call p2p(fldin(:,:,:,4),plin,plout,ps,fldout(:,:,:,3))           !v-wind
  if (trim(qvar) == 'td') then                                           !dew point temperature case
     fldin(:,:,:,5) = calcmet_rh_fromTd(fldin(:,:,:,2),fldin(:,:,:,5))
  else if (trim(qvar) == 'ttd') then                                          !temperature - dew point temperature case
     fldin(:,:,:,5) = calcmet_rh_fromTd(fldin(:,:,:,2),fldin(:,:,:,2)-fldin(:,:,:,5))
  else if (trim(qvar) == 'qq') then                                      !specific humidity case
     fldin(:,:,:,5) = calcmet_rh(fldin(:,:,:,2),fldin(:,:,:,5),plin)
  else if (trim(qvar) /= 'rh') then
     write(*,*) "ERROR!: 'qvar' must be 'rh', 'qq', 'td'. 'qvar' is ", qvar
     call abort
  end if

  call p2p(fldin(:,:,:,5),plin,plout,ps,fldout(:,:,:,4))       
  where(fldout(:,:,:,4)>100.) fldout(:,:,:,4) = 100.
  where(fldout(:,:,:,4)<0.) fldout(:,:,:,4) = 0.
  fldout(:,:,:,4) = calcmet_q(fldout(:,:,:,1), fldout(:,:,:,4), plout(:,:,:))
  if (o3clim) then                                                  !climatorogical ozone
     call p2p(bufr(:,:,:,1),bufr(:,:,:,2),plout,ps,fldout(:,:,:,5),constant_value=0.)
  else                                                              !input ozone
     call p2p(fldin(:,:,:,6),plin,plout,ps,fldout(:,:,:,5),constant_value=0.)         
     fldout(:,:,:,5)=max(0.,fldout(:,:,:,5))
  end if
  call p2p(fldin(:,:,:,7),plin,plout,ps,fldout(:,:,:,6))            !cloud lw
  fldout(:,:,:,4:6)=max(0.,fldout(:,:,:,4:6))
  write(*,*) 'PASS: vertical interpolation'
  write(*,*)

  !------------------------------------
  !create new gfs sigma file header
  !
  write(*,*) 'make new sigma file header ...'
  call make_sighead(head,yr,mn,dy,hr,fhour,jcap,lonr,latr,levs, &
       &             nvcoord,ntrac,idsl,idvc,idvt,idvm,vcoord,idusr=idusr)
  call sigio_aldata(head,data,iret)
  write(*,*) 'PASS: make new sigma file header'
  write(*,*)
  
  !----------------------------------------------
  !harmonic transformation
  write(*,*) 'spherical harmonic transform ...'
  call sptez(0,jcap,4,lonr,latr,data%hs,hs,-1)
  ps = log(ps*1.e-3)  !log(Pa->kPa)
  call sptez(0,jcap,4,lonr,latr,data%ps,ps,-1)
  call sigio_cnvtdv(lonr*latr,lonr*latr,levs,idvc,idvm,ntrac, &
       &            iret,fldout(:,:,:,1),fldout(:,:,:,4:6),head%cpi,-1)
  call sptezm(0,jcap,4,lonr,latr,levs,data%t,fldout(:,:,:,1),-1)
  call sptezmv(0,jcap,4,lonr,latr,levs,data%d,data%z,fldout(:,:,:,2),fldout(:,:,:,3),-1)
  call sptezm(0,jcap,4,lonr,latr,levs,data%q(:,:,1),fldout(:,:,:,4),-1)
  call sptezm(0,jcap,4,lonr,latr,levs,data%q(:,:,2),fldout(:,:,:,5),-1)
  call sptezm(0,jcap,4,lonr,latr,levs,data%q(:,:,3),fldout(:,:,:,6),-1)
  write(*,*) 'PASS: spherical harmonic transform'
  write(*,*)

  !----------------------------
  !write gfs sigma restart file
  !
  write(*,*) 'write sigma file'
  call sigio_rwohdc(luso,'output.siganl',head,data,iret)
  write(*,*) 'PASS:write sigma file'
  write(*,*)

  write(*,*) 'convert complete!'

end program main
