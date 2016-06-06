module flatskyvars
  
  !-----------------------------------------------------------------------
  ! FLATSKY VARIABLES
  !-----------------------------------------------------------------------
  
  ! Map arrays
  real, allocatable :: fsmap(:), fsmapl(:)

  ! Map size and resolutions
  integer npix
  real fov,dpix
  real theta,phi,theta0,phi0,theta1
  real thetamin,thetamax,phimin,phimax
  real thetap,phip,dtheta,dphi

end module flatskyvars
