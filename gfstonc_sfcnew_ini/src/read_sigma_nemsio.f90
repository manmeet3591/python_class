! theia:
! gfortran -O3 -march=native -fPIC -c kinds.f90
! gfortran -O3 -march=native -fPIC -c nemsio_module.f90
! gfortran -O3 -march=native -fPIC -c nemsio_openclose.f90
! gfortran -O3 -march=native -fPIC -c nemsio_read.f90
! f2py -c read_sigma_nemsio.f90 -m read_nemsio --fcompiler=gnu95
! --f90flags='-march=native -O3 -fPIC' kinds.o nemsio_module.o 
! nemsio_openclose.o nemsio_read.o libbacio_4.a libw3nco_d.a
subroutine read_nemsio_header(filename,nlons,nlats,nlevs,idate,nfhour)
   use kinds, only: r_kind
   use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                            nemsio_getheadvar,nemsio_realkind,&
                            nemsio_readrecv,nemsio_init,nemsio_getfilehead
   integer, intent(out) :: nlons,nlats,nlevs,idate(7),nfhour
   character(len=500), intent(in) :: filename
   type(nemsio_gfile) :: gfile
   integer iret
   call nemsio_init(iret=iret)
   if(iret/=0) then
      write(6,*)'problem with nemsio_init, iret=',iret
      stop
   end if
   call nemsio_open(gfile,filename,'READ',iret=iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_open, iret=',iret
      stop
   endif
   call nemsio_getfilehead(gfile,dimx=nlons,dimy=nlats,dimz=nlevs,idate=idate,iret=iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getheadvar, iret=',iret
      stop
   end if
   call nemsio_close(gfile, iret=iret)
end subroutine read_nemsio_header

subroutine read_nemsio_coords(filename, nlons, nlats, nlevs, vcoord, lats, lons)
   use kinds, only: r_kind
   use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                            nemsio_getheadvar,nemsio_realkind,&
                            nemsio_readrecv,nemsio_init,nemsio_getfilehead
   integer, intent(in) :: nlons,nlats,nlevs
   character(len=500), intent(in) :: filename
   type(nemsio_gfile) :: gfile
   real(nemsio_realkind), intent(out) :: vcoord(nlevs+1,3,2)
   real(nemsio_realkind), intent(out) :: lats(nlats),lons(nlons)
   real(nemsio_realkind), dimension(nlons*nlats) :: lons1,lats1
   real(nemsio_realkind), dimension(nlons,nlats) :: lons2,lats2
   integer iret

   call nemsio_init(iret=iret)
   if(iret/=0) then
      write(6,*)'problem with nemsio_init, iret=',iret
      stop
   end if
   call nemsio_open(gfile,filename,'READ',iret=iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_open, iret=',iret
      stop
   endif
   call nemsio_getfilehead(gfile,iret=iret,vcoord=vcoord)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getfilehead (vcoord), iret=',iret
      stop
   endif
   call nemsio_getfilehead(gfile,iret=iret,lat=lats1)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getfilehead (lat), iret=',iret
      stop
   endif
   call nemsio_getfilehead(gfile,iret=iret,lon=lons1)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getfilehead (lon), iret=',iret
      stop
   endif
   call onedtotwod(lats1,lats2,nlons,nlats)
   call onedtotwod(lons1,lons2,nlons,nlats)
   lats = lats2(1,:)
   lons = lons2(:,1)
   call nemsio_close(gfile, iret=iret)
end subroutine read_nemsio_coords

subroutine read_nemsio_griddata(filename, nlons, nlats, nlevs, ug, vg, tempg, zsg, psg, qg, ozg, cwmrg, dpresg, presg)
  use kinds, only: r_kind
  use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                           nemsio_getheadvar,nemsio_realkind,&
                           nemsio_readrecv,nemsio_init,nemsio_getfilehead
  integer, intent(in) :: nlons,nlats,nlevs
  real(r_kind), intent(out),dimension(nlons,nlats,nlevs) :: &
  ug,vg,tempg,qg,ozg,cwmrg,dpresg,presg
  real(r_kind), intent(out), dimension(nlons,nlats) :: psg,zsg
  real(nemsio_realkind) nems_vcoord(nlevs+1,3,2),ak(nlevs+1),bk(nlevs+1)
  real(r_kind), dimension(nlons,nlats) :: press_bot, press_top
  character(len=500), intent(in) :: filename
  real(nemsio_realkind), dimension(nlons*nlats) :: nems_wrk
  real(nemsio_realkind), dimension(nlons,nlats) :: nems_wrk2
  real(r_kind) :: kap,kap1,kapr,rd,cp
  type(nemsio_gfile) :: gfile
  integer iret, idvc
  logical has_dpres
  call nemsio_init(iret=iret)
  if(iret/=0) then
     write(6,*)'problem with nemsio_init, iret=',iret
     stop
  end if
  call nemsio_open(gfile,filename,'READ',iret=iret)
  if (iret/=0) then
     write(6,*)'problem with nemsio_open, iret=',iret
     stop
  endif

  call nemsio_readrecv(gfile,'pres','sfc',1,nems_wrk,iret=iret)
  call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
  psg(:,:) = nems_wrk2
  if (iret/=0) then
      write(6,*)'problem with nemsio_readrecv(ps), iret=',iret
      stop
  endif

  call nemsio_readrecv(gfile,'hgt','sfc',1,nems_wrk,iret=iret)
  call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
  zsg(:,:) = nems_wrk2
  if (iret/=0) then
      write(6,*)'problem with nemsio_readrecv(hgt), iret=',iret
      stop
  endif

  has_dpres = .true.
  call nemsio_readrecv(gfile,'dpres','mid layer',k,nems_wrk,iret=iret)
  if (iret/=0) then
      rd = 287.05; cp = 1004.6
      kap = rd/cp
      kapr = cp/rd
      kap1 = kap + 1.0
      has_dpres = .false.
      call nemsio_getfilehead(gfile,iret=iret,vcoord=nems_vcoord)
      if ( iret /= 0 ) then
         write(6,*) 'problem reading nemsio header ', &
            'vcoord, Status = ',iret
         stop
      endif
      call nemsio_getfilehead(gfile,iret=iret,idvc=idvc)
      if ( idvc == 0 ) then                         ! sigma coordinate, old file format.
         ak = zero
         bk = nems_vcoord(1:nlevs+1,1,1)
      elseif ( idvc == 1 ) then                     ! sigma coordinate
         ak = zero
         bk = nems_vcoord(1:nlevs+1,2,1)
      elseif ( idvc == 2 .or. idvc == 3 ) then      ! hybrid coordinate
         ak = 0.01_r_kind*nems_vcoord(1:nlevs+1,1,1) ! convert to mb
         bk = nems_vcoord(1:nlevs+1,2,1)
      else
         write(6,*)' ***ERROR*** INVALID value for idvc=',idvc
         stop
      endif
      ! pressure thicknesses
      press_bot = psg 
      do k=2,nlevs+1
         press_top = ak(k)+bk(k)*psg
         dpresg(:,:,k-1) = press_bot - press_top
         !print *,k-1,minval(dpresg(:,:,k-1)),maxval(dpresg(:,:,k-1))
         presg(:,:,k-1) = ((press_bot**kap1-press_top**kap1)/&
                          (kap1*(press_bot-press_top)))**kapr
         !print *,k-1,minval(presg(:,:,k-1)),maxval(presg(:,:,k-1))
         press_bot = press_top
      enddo
  endif

  do k=1,nlevs
     if (has_dpres) then
         call nemsio_readrecv(gfile,'dpres','mid layer',k,nems_wrk,iret=iret)
         if (iret /= 0) then
            write(6,*)'problem with nemsio_readrecv(dpress), iret=',iret
            stop
         endif
         call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
         dpresg(:,:,k) = nems_wrk2
         call nemsio_readrecv(gfile,'pres','mid layer',k,nems_wrk,iret=iret)
         if (iret/=0) then
             write(6,*)'problem with nemsio_readrecv(pres), iret=',iret
             stop
         endif
         call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
         presg(:,:,k) = nems_wrk2
     endif
     call nemsio_readrecv(gfile,'ugrd','mid layer',k,nems_wrk,iret=iret)
     if (iret/=0) then
         write(6,*)'problem with nemsio_readrecv(ugrd), iret=',iret
         stop
     endif
     call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
     ug(:,:,k) = nems_wrk2
     call nemsio_readrecv(gfile,'vgrd','mid layer',k,nems_wrk,iret=iret)
     if (iret/=0) then
         write(6,*)'problem with nemsio_readrecv(vgrd), iret=',iret
         stop
     endif
     call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
     vg(:,:,k) = nems_wrk2
     call nemsio_readrecv(gfile,'tmp','mid layer',k,nems_wrk,iret=iret)
     if (iret/=0) then
        write(6,*)'problem with nemsio_readrecv(tmp), iret=',iret
        stop
     endif
     call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
     tempg(:,:,k) = nems_wrk2
     call nemsio_readrecv(gfile,'spfh','mid layer',k,nems_wrk,iret=iret)
     if (iret/=0) then
        write(6,*)'problem with nemsio_readrecv(spfh), iret=',iret
        stop
     endif
     call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
     qg(:,:,k) = nems_wrk2
     call nemsio_readrecv(gfile,'o3mr','mid layer',k,nems_wrk,iret=iret)
     if (iret/=0) then
        write(6,*)'problem with nemsio_readrecv(o3mr), iret=',iret
        stop
     endif
     call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
     ozg(:,:,k) = nems_wrk2
     call nemsio_readrecv(gfile,'clwmr','mid layer',k,nems_wrk,iret=iret)
     if (iret/=0) then
        write(6,*)'problem with nemsio_readrecv(clwmr), iret=',iret
        stop
     endif
     call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
     cwmrg(:,:,k) = nems_wrk2
  enddo
  call nemsio_close(gfile, iret=iret)
end subroutine read_nemsio_griddata

subroutine onedtotwod(data1,data2,nlons,nlats)
   use nemsio_module, only: nemsio_realkind
   integer, intent(in) :: nlons,nlats
   real(nemsio_realkind), intent(in) :: data1(nlons*nlats)
   real(nemsio_realkind), intent(out) :: data2(nlons,nlats)
   integer i,j,n
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      data2(i,j) = data1(n)
   enddo
end subroutine onedtotwod
