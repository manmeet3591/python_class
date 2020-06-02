! theia:
! gfortran -O3 -march=native -fPIC -c kinds.f90
! gfortran -O3 -march=native -fPIC -c nemsio_module.f90
! gfortran -O3 -march=native -fPIC -c nemsio_openclose.f90
! gfortran -O3 -march=native -fPIC -c nemsio_read.f90
! f2py -c read_sfcflx_nemsio.f90 -m read_nemsio --fcompiler=gnu95
! --f90flags='-march=native -O3 -fPIC' kinds.o nemsio_module.o 
! nemsio_openclose.o nemsio_read.o libbacio_4.a libw3nco_d.a
subroutine read_nemsio_header(filename,nlons,nlats,nrecs,idate,nfhour)
   use kinds, only: r_kind
   use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                            nemsio_getheadvar,nemsio_realkind,&
                            nemsio_readrecv,nemsio_init,nemsio_getfilehead
   implicit none
   integer, intent(out) :: nlons,nlats,nrecs,idate(7),nfhour
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
   call nemsio_getheadvar(gfile,'dimx',nlons,iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getheadvar, iret=',iret
      stop
   end if
   call nemsio_getheadvar(gfile,'dimy',nlats,iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getheadvar, iret=',iret
      stop
   end if
   call nemsio_getheadvar(gfile,'nrec',nrecs,iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getheadvar, iret=',iret
      stop
   end if
   call nemsio_getheadvar(gfile,'idate',idate,iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getheadvar, iret=',iret
      stop
   end if
   call nemsio_getheadvar(gfile,'nfhour',nfhour,iret)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getheadvar, iret=',iret
      stop
   end if
   call nemsio_close(gfile, iret=iret)
end subroutine read_nemsio_header

subroutine read_nemsio_latlons(filename, nlons, nlats, lats, lons)
   use kinds, only: r_kind
   use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                            nemsio_getheadvar,nemsio_realkind,&
                            nemsio_readrecv,nemsio_init,nemsio_getfilehead
   implicit none
   integer, intent(in) :: nlons,nlats
   character(len=500), intent(in) :: filename
   type(nemsio_gfile) :: gfile
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
end subroutine read_nemsio_latlons

subroutine read_nemsio_varnames(filename, nrecs, irecnames, ireclevtypes, ireclevs)
   use kinds, only: r_kind
   use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                            nemsio_getheadvar,nemsio_realkind,&
                            nemsio_readrecv,nemsio_init,nemsio_getfilehead
   implicit none
   integer, parameter :: n_str = 32
   integer, intent(in) :: nrecs
   character(len=500), intent(in) :: filename
   character(len=n_str) :: recnames(nrecs),reclevtypes(nrecs)
   integer, intent(out), dimension(nrecs,n_str) :: irecnames
   integer, intent(out), dimension(nrecs,n_str) :: ireclevtypes
   integer, intent(out), dimension(nrecs) :: ireclevs
   type(nemsio_gfile) :: gfile
   integer iret,nrec

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
   call nemsio_getfilehead(gfile,iret=iret,recname=recnames,reclevtyp=reclevtypes, reclev=ireclevs)
   if (iret/=0) then
      write(6,*)'problem with nemsio_getfilehead (recname,reclevtyp), iret=',iret
      stop
   endif
   do nrec=1,nrecs
      call strtoarr(recnames(nrec), irecnames(nrec,:), n_str)
      call strtoarr(reclevtypes(nrec), ireclevtypes(nrec,:), n_str)
   enddo
   call nemsio_close(gfile, iret=iret)
end subroutine read_nemsio_varnames

subroutine read_nemsio_2dgriddata(filename, nlons, nlats, nrecs, irecnames, ireclevtypes,reclevs,grids)
  use kinds, only: r_kind
  use nemsio_module, only: nemsio_gfile,nemsio_open,nemsio_close,&
                           nemsio_getheadvar,nemsio_realkind,&
                           nemsio_readrecv,nemsio_init
  implicit none
  integer, parameter :: n_str = 32
  integer, intent(in) :: nlons,nlats,nrecs
  real(r_kind), intent(out),dimension(nlons,nlats,nrecs) :: grids
  character(len=500), intent(in) :: filename
  integer,  intent(in) :: irecnames(nrecs,n_str),ireclevtypes(nrecs,n_str)
  character(len=n_str) :: recnames(nrecs),reclevtypes(nrecs)
  integer,  intent(in) :: reclevs(nrecs) 
  real(nemsio_realkind), dimension(nlons*nlats) :: nems_wrk
  real(nemsio_realkind), dimension(nlons,nlats) :: nems_wrk2
  type(nemsio_gfile) :: gfile
  integer iret,nrec
  call nemsio_init(iret=iret)
  print *,  'in read_sfcflx_nemsio.f90'
  if(iret/=0) then
     write(6,*)'problem with nemsio_init, iret=',iret
     stop
  end if
  call nemsio_open(gfile,filename,'READ',iret=iret)
  if (iret/=0) then
     write(6,*)'problem with nemsio_open, iret=',iret
     stop
  endif

  do nrec=1,nrecs
      call arrtostr(irecnames(nrec,:),recnames(nrec),n_str)
      call arrtostr(ireclevtypes(nrec,:),reclevtypes(nrec),n_str)
      !print *,nrec,trim(recnames(nrec)),' ',trim(reclevtypes(nrec)), ' ', reclevs(nrec)
      call nemsio_readrecv(gfile,trim(recnames(nrec)),trim(reclevtypes(nrec)),reclevs(nrec),nems_wrk,iret=iret)
      call onedtotwod(nems_wrk,nems_wrk2,nlons,nlats)
      grids(:,:,nrec) = nems_wrk2
      if (iret/=0) then
          write(6,*)'problem with nemsio_readrecv for',trim(recnames(nrec)),' iret=',iret
          stop
      endif
  enddo
end subroutine read_nemsio_2dgriddata

subroutine onedtotwod(data1,data2,nlons,nlats)
   use nemsio_module, only: nemsio_realkind
   implicit none
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
subroutine strtoarr(strin, chararr, n_str)
  integer, intent(in) :: n_str
  character(len=n_str), intent(in) :: strin
  integer, intent(out) ::  chararr(n_str)
  chararr = 32
  do j=1,len_trim(trim(adjustl(strin)))
     chararr(j) = ichar(strin(j:j))
  enddo
end subroutine strtoarr
subroutine arrtostr(chararr, strout, n_str)
  integer, intent(in) :: n_str
  character(len=n_str), intent(out) :: strout
  integer, intent(in) ::  chararr(n_str)
  do j=1,n_str
     strout(j:j) = char(chararr(j))
  enddo
end subroutine arrtostr
