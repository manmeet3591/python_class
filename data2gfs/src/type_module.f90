module type_module

  implicit none

! Symbolic names for kind types of 4-, 2-, and 1-byte integers:

  integer, parameter, public :: i8b = selected_int_kind(15)
  integer, parameter, public :: i4b = selected_int_kind(9)
  integer, parameter, public :: i2b = selected_int_kind(4)

! Symbolic names for kind types of single- and double-precision reals:

  integer, parameter, public :: sp = kind(1.0)
  integer, parameter, public :: dp = kind(1.0d0)

end module type_module

