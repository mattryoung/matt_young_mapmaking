MODULE intreal_types
  ! This module sets the types used in the Fortran 90 modules

  INTEGER, PARAMETER :: i4b = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: sp  = KIND(1.0E0)
  INTEGER, PARAMETER :: dp  = KIND(1.0D0)
  INTEGER, PARAMETER :: lgt = KIND(.TRUE.)
  INTEGER, PARAMETER :: spc = KIND((1.0,  1.0))
  INTEGER, PARAMETER :: dpc = KIND((1.0d0,1.0d0))

  INTEGER,  PARAMETER:: max_i4b = HUGE(1_i4b)
  REAL(sp), PARAMETER:: max_sp  = HUGE(1.0_sp)
  REAL(dp), PARAMETER:: max_dp  = HUGE(1.0_dp)


END MODULE intreal_types
