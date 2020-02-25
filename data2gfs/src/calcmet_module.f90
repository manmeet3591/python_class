module calcmet_module
  use type_module
  use constant_module, only : grav=>earth_gravity, air_cp, gascon=>air_rd
  implicit none
  private

  real(kind=sp), private, parameter :: gamma = 0.0065 

  public :: calcmet_q, calcmet_rh, calcmet_ps, calcmet_z, calcmet_rh_fromTd

  interface calcmet_q
     module procedure calcmet_q_0d, calcmet_q_1d0d, calcmet_q_2d0d, &
          &           calcmet_q_3d1d, calcmet_q_3d3d
  end interface calcmet_q
  interface calcmet_rh
     module procedure calcmet_rh_0d, calcmet_rh_1d0d, calcmet_rh_2d0d, &
          &           calcmet_rh_3d1d, calcmet_rh_3d3d
  end interface calcmet_rh
  interface calcmet_z
     module procedure calcmet_z_3d1d,calcmet_z_3d3d
  end interface calcmet_z
  interface calcmet_rh_fromTd
     module procedure calcmet_rh_fromTd_0d, calcmet_rh_fromTd_3d
  end interface calcmet_rh_fromTd
  
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! function calcmet_q 
  !!  
  !!  calculate specific humidity from relative humidity
  !!
  !!  Argument
  !!    T    :  Temeprature       [K]
  !!    rh   :  relative humidity [%]
  !!    plev :  pressure          [Pa]
  !!
  !!  Return
  !!    q    :  specific humidity [Kg/Kg]
  !!
  !!
  function calcmet_q_0d(T, rh, plev) result (q)    
    real(kind=sp), intent(in) :: T, rh
    real(kind=sp), intent(in) :: plev       !Pa
    real(kind=sp) :: q

    real(kind=sp), parameter :: eps = 0.622 ! Rd/Rv
    real(kind=sp) :: es, e

    es = exp(19.482 - 4303.4 / (T-29.65)) ! in hPa JMA/WMO
    e = es*rh ! x 100 (hPa->Pa) / 100. (%->no dimension)
    q = eps*e / (plev - (1.0-eps)*e)

  end function calcmet_q_0d
  function calcmet_q_1d0d(T, rh, plev) result (q)
    real(kind=sp), dimension(:), intent(in) :: T, rh
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1)) :: q

    integer(kind=i4b) :: i

    do i=1, size(T,1)
       q(i) = calcmet_q_0d(T(i), rh(i), plev)
    end do

  end function calcmet_q_1d0d
  function calcmet_q_2d0d(T, rh, plev) result (q)
    real(kind=sp), dimension(:,:), intent(in) :: T, rh
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2)) :: q

    integer(kind=i4b) :: i, j

    do j=1, size(T,2)
      do i=1, size(T,1)
         q(i,j) = calcmet_q_0d(T(i,j), rh(i,j), plev)
      end do
    end do

  end function calcmet_q_2d0d
  function calcmet_q_3d1d(T, rh, plev) result (q)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, rh
    real(kind=sp), dimension(:), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: q

    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j=1, size(T,2)
          do i=1, size(T,1)
             q(i,j,k) = calcmet_q_0d(T(i,j,k), rh(i,j,k), plev(k))
          end do
       end do
    end do

  end function calcmet_q_3d1d
  function calcmet_q_3d3d(T, rh, plev) result (q)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, rh, plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: q

    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j=1, size(T,2)
          do i=1, size(T,1)
             q(i,j,k) = calcmet_q_0d(T(i,j,k), rh(i,j,k), plev(i,j,k))
             if( q(i,j,k)<0 ) then
             end if
          end do
       end do
    end do

  end function calcmet_q_3d3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! function calcmet_rh
  !!
  !!  calculate relative humidity from T & q
  !!
  !!  Argument
  !!    T    : temperature [K]
  !!    q    : specific humidity [Kg/Kg]
  !!    plev : pressure [Pa]
  !!
  !!  Return
  !!    rh   : relative humidity [%]
  !!
  !!
  function calcmet_rh_0d(T, q, plev) result (rh)
    real(kind=sp), intent(in) :: T, q
    real(kind=sp), intent(in) :: plev
    real(kind=sp) :: rh

    real(kind=sp), parameter :: eps = 0.622 ! Rd/Rv
    real(kind=sp) :: es, e

    es = exp(19.482 - 4303.4 / (T-29.65)) ! in hPa JMA/WMO
    e = q*plev / (eps + (1-eps)*q)
    rh = e/es ! x 100 (hPa->Pa) / 100. (%->no dimension)

  end function calcmet_rh_0d
  function calcmet_rh_1d0d(T, q, plev) result (rh)
    real(kind=sp), dimension(:), intent(in) :: T, q
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1)) :: rh

    integer(kind=i4b) :: i

    do i = 1, size(T,1)
       rh(i) = calcmet_rh_0d(T(i), q(i), plev)
    end do

  end function calcmet_rh_1d0d
  function calcmet_rh_2d0d(T, q, plev) result (rh)
    real(kind=sp), dimension(:,:), intent(in) :: T, q
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2)) :: rh

    integer(kind=i4b) :: i, j

    do j = 1, size(T,2)
       do i = 1, size(T,1)
          rh(i,j) = calcmet_rh_0d(T(i,j), q(i,j), plev)
       end do
    end do

  end function calcmet_rh_2d0d
  function calcmet_rh_3d1d(T, q, plev) result (rh)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, q
    real(kind=sp), dimension(:), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: rh
    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j = 1, size(T,2)
          do i = 1, size(T,1)
             rh(i,j,k) = calcmet_rh_0d(T(i,j,k), q(i,j,k), plev(k))
          end do
       end do
    end do

  end function calcmet_rh_3d1d
  function calcmet_rh_3d3d(T, q, plev) result (rh)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, q
    real(kind=sp), dimension(:,:,:), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: rh
    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j = 1, size(T,2)
          do i = 1, size(T,1)
             rh(i,j,k) = calcmet_rh_0d(T(i,j,k), q(i,j,k), plev(i,j,k))
          end do
       end do
    end do

  end function calcmet_rh_3d3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! function calcmet_ps
  !!
  !!  calcumate surface pressure
  !!
  !!  Argument
  !!    p  : pressure [Pa]
  !!    T  : temperature [T]
  !!    z  : geopotential height [m]
  !!    zs : geometirc height [m]
  !!
  !!  Return
  !!    ps : surface pressure [Pa]
  !!
  !!
  function calcmet_ps(p, T, z, zs) result(ps)
    real(kind=sp), dimension(:,:,:), intent(in) :: p, T, z
    real(kind=sp), dimension(:,:), intent(in) :: zs
    real(kind=sp), dimension(size(zs,1), size(zs,2)) :: ps

    real(kind=sp), parameter :: ggg = grav/gascon/gamma
    integer(kind=i4b), dimension(1) :: kk
    integer(kind=i4b) :: i, j, k

    do j=1, size(zs,2)
      do i=1, size(zs, 1)
        kk = minloc(abs(zs(i,j)-z(i,j,:)),mask=zs(i,j)<z(i,j,:))
        k = kk(1)
        ps(i,j) = p(i,j,k) * (1.0 - gamma*(zs(i,j)-z(i,j,k))/T(i,j,k))**ggg
      end do
    end do
  end function calcmet_ps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! function calcmet_z
  !!
  !!  calculate geopotential height
  !!
  !!  Argument
  !!    T  : temperature [K]
  !!    zs : geometirc height [m]
  !!    ps : surface pressure [Pa]
  !!    pl : pressure [Pa]
  !!
  !!  Return
  !!    z  : geopotential height [m]
  !!
  function calcmet_z_3d1d(T,zs,ps,pl) result(z)
    real(kind=sp), dimension(:,:,:), intent(in) :: T !temperature
    real(kind=sp), dimension(:,:), intent(in) :: zs,ps  !surface height, surface pressure
    real(kind=sp), dimension(:), intent(in) :: pl    !pressure 

    real(kind=sp), parameter :: gg = gascon/grav
    real(kind=sp), dimension(size(T,1),size(T,2),size(T,3)) :: z

    integer(kind=i4b) :: k, n

    n = size(pl)
    z(:,:,1) = zs(:,:) - gg*T(:,:,1)*log(pl(1)/ps(:,:))
!print *, "zs=",zs(240,120), " T=", T(240,120,1), " lns(1)=", lns(1), " z=",z(240,120,1)
    do k=2, n
      z(:,:,k) = z(:,:,k-1) + 0.5*gg*(T(:,:,k-1)+T(:,:,k))*(log(pl(k-1))-log(pl(k)))
    end do

  end function calcmet_z_3d1d
  function calcmet_z_3d3d(T,zs,ps,pl) result(z)
    real(kind=sp), dimension(:,:,:), intent(in) :: T  !temperature
    real(kind=sp), dimension(:,:), intent(in) :: zs,ps   !surface height, surface pressure
    real(kind=sp), dimension(:,:,:), intent(in) :: pl !pressure 

    real(kind=sp), parameter :: gg = gascon/grav
    real(kind=sp), dimension(size(T,1),size(T,2),size(T,3)) :: z

    integer(kind=i4b) :: k, n

    n = size(pl,3)
    z(:,:,1) = zs(:,:) - gg*T(:,:,1)*log(pl(:,:,1)/ps(:,:))
!print *, "zs=",zs(240,120), " T=", T(240,120,1), " lns(1)=", lns(1), " z=",z(240,120,1)
    do k=2, n
      z(:,:,k) = z(:,:,k-1) + 0.5*gg*(T(:,:,k-1)+T(:,:,k))*(log(pl(:,:,k-1))-log(pl(:,:,k)))
    end do

  end function calcmet_z_3d3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! function calcmet_rh_fromTd
  !!
  !!  calculate relative humidity from dew point temperature
  !!
  !!  Argument
  !!    T  : temperature [K]
  !!    Td : dew pointe temperature [K]
  !!
  !!  Return
  !!    rh : relative humidity [%]
  !!
  !!
  function calcmet_rh_fromTd_0d(T, Td) result (rh)
    real(kind=sp), intent(in) :: T, Td
    real(kind=sp) :: rh

    real(kind=sp), parameter :: eps = 0.622 ! Rd/Rv
    real(kind=sp) :: es, e

    es = exp(19.482 - 4303.4 / (T-29.65))  ! in hPa JMA/WMO
    e  = exp(19.482 - 4303.4 / (Td-29.65)) ! in hPa JMA/WMO
    rh = e/es * 100 

  end function calcmet_rh_fromTd_0d
  function calcmet_rh_fromTd_3d(T, Td) result (rh)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, Td
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: rh
    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j = 1, size(T,2)
          do i = 1, size(T,1)
             rh(i,j,k) = calcmet_rh_fromTd_0d(T(i,j,k), Td(i,j,k))
          end do
       end do
    end do

  end function calcmet_rh_fromTd_3d

end module calcmet_module
