program make_cltable

  use cltable
  use cosmology

  implicit none

  ! So far 3 models implemented
  !    model = 1 --> BPPS2 AGN Feedback Delta = 200
  !    model = 2 --> BPPS2 AGN Feedback Delta = 500
  !    model = 3 --> BPPS2    Adiabatic Delta = 500

  integer model

  ! Output filename
  character *128 fileout

  integer i4arg

  ! Usage check
  if(iargc()/=2) then
     write(*,11) 
11   format('Usage: make_cltable <fileout> <model>')
     stop
  endif

  ! Get commandline
  call getarg(1,fileout)
  model = i4arg(2,1)

  ! Set background cosmology
  omegam = 0.25
  omegab = 0.043
  omegal = 0.75
  h = 0.7
  sigma8 = 0.8
  ns = 0.96
  w = -1
  fb = omegab/omegam

  rho_0     = 2.775e11*omegam*0.7**2
  rhocrit_0 = 2.775e11*0.7**2  

  ! Make table 
  call makecltable(fileout,model)

  stop
  
end program make_cltable
