! theia:
! gfortran -O3 -march=native -O3 -fPIC -c shtns.f90
! gfortran -O3 -march=native -O3 -fPIC -c kinds.f90
! gfortran -O3 -march=native -O3 -fPIC -c sigio_module.f90
! f2py -c read_sigma.f90 -m read_sigma --fcompiler=gnu95
! --f90flags='-march=native -O3 -fopenmp' -lgomp shtns.o kinds.o sigio_module.o
! -L/scratch3/BMC/gsienkf/whitaker/contrib/pythonextras/Anaconda-1.8.0/lib
! -lshtns -lfftw3
subroutine get_vcoord(filename, nlevs, nvcoord, vcoord)
  use sigio_module, only: sigio_sclose,sigio_swohdc,sigio_head,sigio_srhead,&
  sigio_srohdc,sigio_aldata,sigio_data,sigio_sropen,sigio_srdata,sigio_axdata
  integer, intent(in) :: nlevs, nvcoord
  real, intent(out) :: vcoord(nlevs+1,nvcoord)
  character(len=500), intent(in) :: filename
  type(sigio_head) sighead
  integer lu,iret
  lu = 7
  call sigio_sropen(lu,trim(filename),iret)
  if (iret .ne. 0) then
     print *,'error opening ',trim(filename),iret
     stop
  endif
  call sigio_srhead(lu,sighead,iret)
  if (iret .ne. 0) then
     print *,'error reading header from ',trim(filename),iret
     stop
  else
     vcoord = sighead%vcoord(:,:)
  endif 
  call sigio_sclose(lu,iret)
end subroutine get_vcoord
subroutine read_header(filename, nlons, nlats, nlevs, ntrunc, nvcoord, idate,&
  fhour)
  use sigio_module, only: sigio_sclose,sigio_swohdc,sigio_head,sigio_srhead,&
  sigio_srohdc,sigio_aldata,sigio_data,sigio_sropen,sigio_srdata,sigio_axdata
  implicit none
  integer, intent(out) :: nlons,nlats,nlevs,ntrunc,nvcoord,idate(4)
  real, intent(out) :: fhour
  character(len=500), intent(in) :: filename
  type(sigio_head) sighead
  integer lu,iret
  lu = 7
  call sigio_sropen(lu,trim(filename),iret)
  if (iret .ne. 0) then
     print *,'error opening ',trim(filename),iret
     stop
  endif
  call sigio_srhead(lu,sighead,iret)
  if (iret .ne. 0) then
     print *,'error reading header from ',trim(filename),iret
     stop
  else
     nlons = sighead%lonb
     nlats = sighead%latb
     nlevs = sighead%levs
     ntrunc = sighead%jcap
     nvcoord = sighead%nvcoord
     idate = sighead%idate
     fhour = sighead%fhour
  endif 
  call sigio_sclose(lu,iret)
end subroutine read_header
subroutine read_specdata(filename, ntrunc, nlevs, vrtspec, divspec, virtempspec, topospec, lnpsspec, spfhumspec,&
  ozspec,cwmrspec)
  use sigio_module, only: sigio_sclose,sigio_swohdc,sigio_head,sigio_srhead,&
  sigio_srohdc,sigio_aldata,sigio_data,sigio_sropen,sigio_srdata,sigio_axdata
  use kinds, only: r_kind
  implicit none
  complex(r_kind), intent(out),dimension((ntrunc+1)*(ntrunc+2)/2,nlevs) :: &
  vrtspec, divspec, virtempspec, spfhumspec, ozspec, cwmrspec
  complex(r_kind), intent(out),dimension((ntrunc+1)*(ntrunc+2)/2) :: &
  topospec, lnpsspec
  integer, intent(in) :: ntrunc,nlevs
  real(r_kind) norm
  character(len=500), intent(in) :: filename
  integer lu,k,iret,ndimspec
  type(sigio_data) :: sigdata
  type(sigio_head) :: sighead
  lu = 7
  !print *,trim(filename)
  call sigio_srohdc(lu,trim(filename),sighead,sigdata,iret)
  if (iret .ne. 0) then
    print *,'error reading ',trim(filename),iret
    stop
  endif
  ndimspec = (ntrunc+1)*(ntrunc+2)/2
  norm=1
  call copyspecin(sigdata%ps, lnpsspec,ndimspec,norm)
  call copyspecin(sigdata%hs, topospec,ndimspec,norm)
  do k=1,nlevs
     call copyspecin(sigdata%z(:,k),vrtspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%d(:,k),divspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%t(:,k),virtempspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%q(:,k,1),spfhumspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%q(:,k,2),ozspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%q(:,k,3),cwmrspec(:,k),ndimspec,norm)
  enddo
  call sigio_axdata(sigdata,iret)
  call sigio_sclose(lu,iret)
end subroutine read_specdata
subroutine read_griddata(filename, nlons, nlats, nlevs, ug, vg, tempg, zsg, psg, qg, ozg, cwmrg)
  use shtns, only: getuv, grdtospec, spectogrd, shtns_init, nlm
  use sigio_module, only: sigio_sclose,sigio_swohdc,sigio_head,sigio_srhead,&
  sigio_srohdc,sigio_aldata,sigio_data,sigio_sropen,sigio_srdata,sigio_axdata
  use kinds, only: r_kind
  implicit none
  real(r_kind),parameter:: rerth  =6.3712e+6      ! radius of earth   (m)
  real(r_kind), intent(out),dimension(nlons,nlats,nlevs) :: &
  ug,vg,tempg,qg,ozg,cwmrg
  real(r_kind), intent(out), dimension(nlons,nlats) :: psg,zsg
  integer, intent(in) :: nlons,nlats,nlevs
  character(len=500), intent(in) :: filename
  real(r_kind) norm
  integer lu,k,iret,ntrunc,ndimspec
  complex(r_kind), allocatable, dimension(:,:) :: vrtspec,divspec,virtempspec,spfhumspec,&
  ozspec,cwmrspec
  complex(r_kind), allocatable, dimension(:) :: topospec,lnpsspec
  type(sigio_data) :: sigdata
  type(sigio_head) :: sighead
  lu = 7
  !print *,trim(filename)
  call sigio_srohdc(lu,trim(filename),sighead,sigdata,iret)
  if (iret .ne. 0) then
    print *,'error reading ',trim(filename),iret
    stop
  endif
  ntrunc = sighead%jcap
  ndimspec = (ntrunc+1)*(ntrunc+2)/2
  ! initialize spherical harmonic lib
  call shtns_init(nlons,nlats,ntrunc)
  !print *,'nlm,ndimspec',nlm,ndimspec
  ! convert spectral arrays to double precision complex,
  ! re-normalize coefficients.
  allocate(vrtspec(ndimspec,nlevs))
  allocate(divspec(ndimspec,nlevs))
  allocate(virtempspec(ndimspec,nlevs))
  allocate(spfhumspec(ndimspec,nlevs))
  allocate(ozspec(ndimspec,nlevs))
  allocate(cwmrspec(ndimspec,nlevs))
  allocate(lnpsspec(ndimspec),topospec(ndimspec))
  norm = sqrt(2.*4.*atan(1.0))
  call copyspecin(sigdata%ps, lnpsspec,ndimspec,norm)
  call copyspecin(sigdata%hs, topospec,ndimspec,norm)
!$omp parallel do private(k)
  do k=1,nlevs
     call copyspecin(sigdata%z(:,k),vrtspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%d(:,k),divspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%t(:,k),virtempspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%q(:,k,1),spfhumspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%q(:,k,2),ozspec(:,k),ndimspec,norm)
     call copyspecin(sigdata%q(:,k,3),cwmrspec(:,k),ndimspec,norm)
  enddo
!$omp end parallel do
  call spectogrd(lnpsspec, psg)
  call spectogrd(topospec, zsg)
  if ( mod(sighead%idvm,10) /= 2 ) then
    psg = 1000.*exp(psg) ! convert to Pa
  else
    psg = 1000.*psg ! convert to Pa
  endif
  !print *,'min/max ps',minval(psg),maxval(psg)
  !print *,'min/max zs',minval(zsg),maxval(zsg)
!$omp parallel do private(k)
  do k=1,nlevs
     call getuv(vrtspec(:,k),divspec(:,k),ug(:,:,k),vg(:,:,k),rerth)
     call spectogrd(virtempspec(:,k),tempg(:,:,k))
     call spectogrd(spfhumspec(:,k),qg(:,:,k))
     call spectogrd(ozspec(:,k),ozg(:,:,k))
     call spectogrd(cwmrspec(:,k),cwmrg(:,:,k))
  enddo
!$omp end parallel do
  call sigio_axdata(sigdata,iret)
  call sigio_sclose(lu,iret)
  deallocate(lnpsspec,topospec)
  deallocate(vrtspec,divspec,virtempspec,spfhumspec,ozspec,cwmrspec)
end subroutine read_griddata
subroutine copyspecin(rspecdata,cspecdata,ndimspec,norm)
  use kinds, only: r_kind, r_single
  implicit none
  real(r_kind), intent(in) :: norm
  integer, intent(in) :: ndimspec
  real(r_single), intent(in) :: rspecdata(2*ndimspec)
  complex(r_kind), intent(out) :: cspecdata(ndimspec)
  integer n,nn
  nn = 1
  ! factor of sqrt(2.*pi) accounts for difference in normalization
  ! between ncep libs and shtns (which uses orthonormalized norm)
  do n=1,ndimspec
     cspecdata(n) = norm*cmplx(rspecdata(nn),rspecdata(nn+1))
     nn = nn + 2
  enddo
end subroutine copyspecin
