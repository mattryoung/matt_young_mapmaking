module array_bounds

  integer, parameter :: n1         =  256
  integer, parameter :: n2         =  n1
  integer, parameter :: n3         =  n1
  integer, parameter :: np1d       =  n1
  integer, parameter :: np2d       =  n1
  integer, parameter :: np3d       =  n1
  integer, parameter :: nFmax      =  n1*n2*n3
  integer, parameter :: ircmax     =  max(n1,n2,n3)/2

  integer, parameter :: itabmax    =  2001, &
#ifdef DARWIN
                        np         =  11000000
  integer, parameter :: nodemax    =  40000000
#else
                        np         =  110000000
  integer, parameter :: nodemax    =  400000000
#endif
  integer, parameter :: Npv_evmax =   10,   &
                        Npv_evmin  = -10,   &
                        Nevmax     =  19,   &
                        Nfscmax    =  10
  integer, parameter :: nbormax    =  500000
  integer, parameter :: ntabmax    =  2000


end module array_bounds
