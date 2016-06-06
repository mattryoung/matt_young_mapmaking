program make_maptable

  use maptable
  use cosmology
  use textlib

  implicit none

  ! So far 3 models implemented
  !    model = 1 --> BPPS2 AGN Feedback Delta = 200
  !    model = 2 --> BPPS2 AGN Feedback Delta = 500
  !    model = 3 --> BPPS2    Adiabatic Delta = 500

  integer model
  real nu, alphaCMB
  ! Output filename
  character *128 fileout

  ! Usage check
  if(iargc()<3) then
     write(*,11) 
11   format('Usage: make_maptable <fileout> <model> [<frequency> <alpha_CMB>]')
     stop
  endif

  ! Get commandline
  call getarg(1,fileout)
  model = i4arg(2,1)
  nu    = r4arg(3,0.0)
  alphaCMB = r4arg(4,0.0)

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
  call makemaptable(fileout,model,nu,alphaCMB)

  stop
  
end program make_maptable
