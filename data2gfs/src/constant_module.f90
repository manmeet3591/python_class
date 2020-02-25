module constant_module
  use type_module
  implicit none
  public

  !math constants
  real(kind=dp), parameter :: math_pi  = 3.14159265358979323846_dp

  !Geophysics constants
  real(kind=dp), parameter :: earth_radius = 6.371e6_dp             !radius of earth(m)
  real(kind=dp), parameter :: earth_daysec = 86400.0_dp             !seconds per day(s)
  real(kind=dp), parameter :: earth_omega  = 2.0_dp*math_pi/earth_daysec !ang vel of earth(1/s)
  real(kind=dp), parameter :: earth_gravity = 9.80665_dp            !gravity (ms**-2)

  !air constants
  real(kind=dp), parameter :: air_rd = 287.0_dp    ! gas constant of dry air, J/deg/kg
  real(kind=dp), parameter :: air_cp = 1004.0_dp   ! specific heat at constant pressure, J/deg/kg
  real(kind=dp), parameter :: air_cv =  717.0_dp   ! specific heat at constant volume, J/deg/kig
  real(kind=dp), parameter :: air_kappa = air_rd/air_cp
  real(kind=dp), parameter :: air_gamma = air_cp/air_cv

end module constant_module
