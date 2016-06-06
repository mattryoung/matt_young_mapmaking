module intreal_types

  use, intrinsic :: iso_c_binding

  !-----------------------------------------------------------------------
  ! TYPE DEFINITIONS
  !-----------------------------------------------------------------------

  integer, parameter :: dp=selected_real_kind(12,200)
  integer, parameter :: sp=selected_real_kind(5,30)
  integer, parameter :: i4b = SELECTED_INT_KIND(9)
  integer, parameter :: i8b=selected_int_kind(16)
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
  integer, parameter :: spc=kind((1.0_sp,1.0_sp))
  integer, parameter :: lgt = KIND(.TRUE.)

  
end module intreal_types
