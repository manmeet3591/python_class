module gfs_module
  use type_module
  use sigio_module
  use sigio_r_module
  use ip_module
  implicit none
  private

  public :: make_sighead

contains
  subroutine make_sighead(head,yr,mn,dy,hr,fhour,jcap,lonr,latr,levs, &
       &                   nvcoord,ntrac,idsl,idvc,idvt,idvm,vcoord,idusr)
    type(sigio_head), intent(inout) :: head
    integer(kind=i4b), intent(in)   :: yr,mn,dy,hr
    real(kind=sp), intent(in)       :: fhour
    integer(kind=i4b), intent(in)   :: jcap,lonr,latr,levs,nvcoord,ntrac,idsl,idvc,idvt,idvm
    real(kind=sp), intent(in)       :: vcoord(:,:)
    integer(kind=sp), intent(in),optional :: idusr
    integer(kind=i4b) :: iret

    head%nvcoord = nvcoord
    head%levs = levs
    head%ntrac = ntrac
    head%idvm = idvm
    call sigio_alhead(head,iret,levs,nvcoord,ntrac,idvm)
    if (nvcoord == 3) then
       head%ivs = 200509
    else
       head%ivs = 198410
    end if
    !  head%clabsig=
    head%fhour = fhour
    head%idate(1) = hr
    head%idate(2) = mn
    head%idate(3) = dy
    head%idate(4) = yr
    head%jcap = jcap
    head%levs = levs
    head%vcoord = vcoord
    head%itrun  = 1  !truncation flag (1=triangular)
    head%iorder = 2  !coefficient order flag (2=IBM order)
    head%irealf = 1  !floating point flag (1=4-byte ieee)
    head%igen = 255  !model generating flag (see ON384) (not used in GFS)
    head%latf = latr
    head%lonf = lonr
    head%latb = latr
    head%lonb = lonr
    head%latr = latr
    head%lonr = lonr
    head%ntrac = ntrac
    head%icen2 = 0    !???????
    head%iens = (0,0) !???????
    head%idpp = 0     !???????
    head%idsl = idsl    !semi-lagrangian id (1=Phillips 2=middle of layer)
    head%idvc = idvc    !vertical coordinate (1=sigma 2=hybrid 3=general hybrid)
    head%idvm = idvm    !thermodynamic id
    head%idvt = idvt   !tracer variable id (21=vapor, o3, cloud)
    head%idrun = 0   !run id
    head%idusr = 0   !user defiend id
    head%ncldt = 1   !number of cloud type
    !head%cfvars(:)
    head%nxgr = 0
    head%nxss = 0
    if (present(idusr)) then
       head%idusr = idusr   !user defiend id
    end if

  end subroutine make_sighead

end module gfs_module
