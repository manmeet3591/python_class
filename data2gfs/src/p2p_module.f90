module p2p_module
  use type_module
  use constant_module, only : g=>earth_gravity, air_rd, air_cp
  implicit none
  private

  public :: p2p, p2p_extrapolate_T
  
  interface idxsearch
     module procedure idxsearch_0d, idxsearch_1d
  end interface idxsearch

contains
   subroutine p2p(in,pin,pout,ps,out,constant_value)
    !!
    !! vertical interpolation from input pressure level to output pressure level 
    !! by log-p cubic interpolation method
    !!
    !! CAUTION
    !!  - the value under the lowest input level or the surface level is  
    !!    extrapolated by constant value of the lowest input level or
    !!    'constant_value' if specified.
    !!    
    real(kind=sp), intent(in)           :: in(:,:,:)      !input fieled
    real(kind=sp), intent(in)           :: pin(:,:,:)     !input pressure filed
    real(kind=sp), intent(in)           :: pout(:,:,:)    !output pressure field
    real(kind=sp), intent(in)           :: ps(:,:)        !surface pressure
    real(kind=sp), intent(out)          :: out(:,:,:)     !output fieled
    real(kind=sp), intent(in), optional :: constant_value !constant value under lowest input value

    integer(kind=i4b) :: i,j,k,ko,xn,yn,k1,k2,kl
    integer(kind=i4b) :: idx(size(out,3)) 
    real(kind=sp) :: po,p1,p2,p3,p4,q1,q2,q3,q4,a
    real(kind=sp) :: lnpi(size(in,3)),lnpo(size(out,3))

    xn = size(in,1)
    yn = size(in,2)
    k1 = size(in,3)
    k2 = size(out,3)

    !$omp parallel private(j,k,kl,ko,lnpi,lnpo,idx,a,po,p1,p2,p3,p4,q1,q2,q3,q4)
    !$omp do
    do i = 1, xn
       do j = 1, yn
          lnpi = -log(pin(i,j,:))
          call idxsearch(pin(i,j,:), ps(i,j), ko) !search for surface pressure level index
          kl = ko + 1
          lnpo = -log(pout(i,j,:))
          call idxsearch(lnpi, lnpo, idx)
          do k = 1, k2
             ko = idx(k)            
             if (ko == 0 .or. pin(i,j,ko) > ps(i,j)) then !outside of lowest level
                if (present(constant_value)) then
                   out(i,j,k) = constant_value
                else
                   out(i,j,k) = in(i,j,ko+1)
                end if
             else if (ko == k1) then !outside top level is constant
                out(i,j,k) = in(i,j,k1)
             else if (ko == kl .or. ko+1 == k1) then !linear interpolatlion
                a = (in(i,j,ko+1) - in(i,j,ko))/(lnpi(ko+1) - lnpi(ko))
                out(i,j,k) = a*(lnpo(k) - lnpi(ko)) + in(i,j,ko)
             else                    !cubic lagrangian x1<x2<xo<x3<x4
                po = lnpo(k)
                p1 = lnpi(ko-1)
                p2 = lnpi(ko)
                p3 = lnpi(ko+1)
                p4 = lnpi(ko+2)
                
                q1 = (po-p2)*(po-p3)*(po-p4)/(p1-p2)/(p1-p3)/(p1-p4)
                q2 = (po-p1)*(po-p3)*(po-p4)/(p2-p1)/(p2-p3)/(p2-p4)
                q3 = (po-p1)*(po-p2)*(po-p4)/(p3-p1)/(p3-p2)/(p3-p4)
                q4 = (po-p1)*(po-p2)*(po-p3)/(p4-p1)/(p4-p2)/(p4-p3)
                
                out(i,j,k) = in(i,j,ko-1)*q1 + in(i,j,ko)*q2 +  &
                     &       in(i,j,ko+1)*q3 + in(i,j,ko+2)*q4
             end if
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine p2p

   subroutine p2p_extrapolate_T(tin,pin,pout,psin,psout,zs,tout,constant_value)
    !!
    !! vertical interpolation from input pressure level to output pressure level 
    !! by log-p cubic interpolation method
    !!
    !! CAUTION
    !!  - the value under the lowest input level or the surface level is  
    !!    extrapolated by adiabatic lapse late
    !!    
    real(kind=sp), intent(in)           :: tin(:,:,:)     !input temperature fieled
    real(kind=sp), intent(in)           :: pin(:,:,:)     !input  pressure filed
    real(kind=sp), intent(in)           :: pout(:,:,:)    !output pressure field
    real(kind=sp), intent(in)           :: psin(:,:)      !input surface pressure [Pa]
    real(kind=sp), intent(in)           :: psout(:,:)     !output surface pressure [Pa]
    real(kind=sp), intent(in)           :: zs(:,:)        !output surface height [m]
    real(kind=sp), intent(out)          :: tout(:,:,:)    !output temperature fieled
    real(kind=sp), intent(in), optional :: constant_value !constant value under lowest input value

    integer(kind=i4b) :: i,j,k,ko,xn,yn,k1,k2,kl
    integer(kind=i4b) :: idx(size(tout,3)) 
    real(kind=sp) :: po,p1,p2,p3,p4,q1,q2,q3,q4,a,Ts
    real(kind=sp) :: lnpi(size(tin,3)),lnpo(size(tout,3))

    xn = size(tin,1)
    yn = size(tin,2)
    k1 = size(tin,3)
    k2 = size(tout,3)

    !$omp parallel private(j,k,kl,ko,lnpi,lnpo,idx,a,po,p1,p2,p3,p4,q1,q2,q3,q4,Ts)
    !$omp do
    do i = 1, xn
       do j = 1, yn
          lnpi = -log(pin(i,j,:))
          call idxsearch(pin(i,j,:), psin(i,j), ko)
          kl = ko + 1
          Ts = extrapolate_Ts(tin(i,j,kl), pin(i,j,kl), psout(i,j))
          lnpo = -log(pout(i,j,:))
          call idxsearch(lnpi, lnpo, idx)
          do k = 1, k2
             ko = idx(k)            
             if (ko == 0 .or. pin(i,j,ko) > psin(i,j)) then !outside of lowest level
                if (present(constant_value)) then
                   tout(i,j,k) = constant_value
                else
                   tout(i,j,k) = extrapolate_T(pout(i,j,k),psout(i,j),zs(i,j),Ts)
                end if
             else if (ko == k1) then !outside top level is constant
                tout(i,j,k) = tin(i,j,k1)
             else if (ko == kl .or. ko+1 == k1) then !linear interpolatlion
                a = (tin(i,j,ko+1) - tin(i,j,ko))/(lnpi(ko+1) - lnpi(ko))
                tout(i,j,k) = a*(lnpo(k) - lnpi(ko)) + tin(i,j,ko)
             else                    !cubic lagrangian x1<x2<xo<x3<x4
                po = lnpo(k)
                p1 = lnpi(ko-1)
                p2 = lnpi(ko)
                p3 = lnpi(ko+1)
                p4 = lnpi(ko+2)
                
                q1 = (po-p2)*(po-p3)*(po-p4)/(p1-p2)/(p1-p3)/(p1-p4)
                q2 = (po-p1)*(po-p3)*(po-p4)/(p2-p1)/(p2-p3)/(p2-p4)
                q3 = (po-p1)*(po-p2)*(po-p4)/(p3-p1)/(p3-p2)/(p3-p4)
                q4 = (po-p1)*(po-p2)*(po-p3)/(p4-p1)/(p4-p2)/(p4-p3)
                
                tout(i,j,k) = tin(i,j,ko-1)*q1 + tin(i,j,ko)*q2 +  &
                     &       tin(i,j,ko+1)*q3 + tin(i,j,ko+2)*q4
             end if
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine p2p_extrapolate_T

  subroutine idxsearch_0d(ar1, val, idx)
    !ar1 assumed to be monotonically decreasing or increasing
    !idx=k means val is between ar1(k) and ar1(k+1)
    !idx=0 or idx=size(ar1) indicates val is outside the range of ar1

    real(kind=sp), intent(in)      :: ar1(:) !value to search
    real(kind=sp), intent(in)      :: val    !set of value to serach for
    integer(kind=i4b), intent(out) :: idx    !interval locations 

    integer(kind=i4b) :: k1, j

    k1 = size(ar1)

    if (ar1(1) < ar1(k1)) then
       !monotonically increasing case
       if (val <= ar1(1)) then
          idx = 0
          return
       else if ( val >= ar1(k1)) then
          idx = k1
          return
       end if
       do j = 1, k1-1
          if (val <= ar1(j+1)) then
             idx = j
             return
          end if
       end do
    else
       !monotonically decreasing case
       if (val >= ar1(1)) then
          idx = 0
          return
       else if ( val <= ar1(k1)) then
          idx = k1
          return
       end if
       do j = 1, k1-1
          if (val >= ar1(j+1)) then
             idx = j
             return
          end if
       end do
    end if

  end subroutine idxsearch_0d
  subroutine idxsearch_1d(ar1, ar2, idx)
    !ar1 and ar2  assumed to be monotonically decreasing or increasing
    !idx(k) means ar2(k) is between ar1(k) and ar1(k+1)
    !idx(k) = 0 or size(ar1) indicates ar2(k) is outside the range of ar1

    real(kind=sp), intent(in)      :: ar1(:) !sequence value to search
    real(kind=sp), intent(in)      :: ar2(:) !set of value to serach for
    integer(kind=i4b), intent(out) :: idx(:) !interval locations 

    integer(kind=i4b) :: k1, k2, i, j
    real(kind=sp) :: val

    k1 = size(ar1)
    k2 = size(ar2)

    if (ar1(1) < ar1(k1)) then
       !monotonically increasing case
       do i = 1, k2
          val = ar2(i)
          if (val <= ar1(1)) then
             idx(i) = 0
             cycle
          else if ( val >= ar1(k1)) then
             idx(i) = k1
             cycle
          end if
          do j = 1, k1-1
             if (val <= ar1(j+1)) then
                idx(i) = j
                exit
             end if
          end do
       end do
    else
       !monotonically decreasing case
       do i = 1, k2
          val = ar2(i)
          if (val >= ar1(1)) then
             idx(i) = 0
             cycle
          else if ( val <= ar1(k1)) then
             idx(i) = k1
             cycle
          end if
          do j = 1, k1-1
             if (val >= ar1(j+1)) then
                idx(i) = j
                exit
             end if
          end do
       end do
    end if

  end subroutine idxsearch_1d

  function extrapolate_T(pl, ps, zs, Ts) result(T)
    real(kind=sp), intent(in) :: pl     !pressure
    real(kind=sp), intent(in) :: ps     !surface pressure
    real(kind=sp), intent(in) :: zs     !surface height
    real(kind=sp), intent(in) :: Ts     !surface temperature

    real(kind=sp), parameter :: zc1 = 2000., zc2 = 2500., Tc = 298., dTdz=0.0065
    real(kind=sp) :: T, gamma, T0, T1, y    

    if (zs<zc1) then
      gamma = dTdz
    else
      T1 = Ts + dTdz*zs
      T0 = min(T1, Tc) ! value for zs>zc2
      if (zs<=zc2) then
        T0 = (T0-T1)/(zc2-zc1)*(zs-zc1) + T1
      end if
      gamma = max(T0-Ts,0.)/zs
    end if
    y = gamma*air_rd/g*log(pl/ps)
    T = Ts*(1 + y + y*y/2 + y**3/6)

  end function extrapolate_T

  function extrapolate_Ts(Tl, pl, ps) result(Ts)
    real(kind=sp), intent(in) :: Tl   !lowest level temperature
    real(kind=sp), intent(in) :: pl   !lowest level pressure
    real(kind=sp), intent(in) :: ps   !surface pressure
    real(kind=sp) :: Ts
    real(kind=sp), parameter :: dTdz=0.0065 !K/m

    Ts = Tl*(1 + dTdz*air_rd/g*(1./pl*ps-1.))

  end function extrapolate_Ts

end module p2p_module
