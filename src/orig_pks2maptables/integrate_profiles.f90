module integrate_profiles

  use cosmology
  use bbps_profile

  implicit none

  real, parameter :: ds = 1e-2, dx = 1e-2, xmax = 4

contains

  real function integrate_ptilde(theta,z,mh)

    implicit none

    integer i,n
    real theta,z,mh
    real Delta,chi,rvir,b,s,r,x,s0
    real y,rmax

    rmax = 4

    integrate_ptilde = 0

    Delta = bbps_delta(bbps_model)

    chi = rofz(z)
    rvir = (3*mh/4/pi/Delta/rhocrit(z))**(1./3.)     

    b = chi * sin(theta) / rvir 

    if(b>rmax) return

!    s0 = sqrt(rmax**2-b**2) ! this is for a spherical cut at rmax
    s0 = rmax                ! this is for a cylindrical cut at rmax

    n = int(s0/ds) + 1

    y = 0
    do i=1,n
       s = (i-0.5) * ds
       x = sqrt(s**2+b**2)
       y = y + bbps_ptilde(mh,z,x)
    enddo

    y = y * ds

    integrate_ptilde = y ! this thing is dimensionless

    return 

  end function integrate_ptilde

  real function integrate_rhotilde(theta,z,mh)

    implicit none

    integer i,n
    real theta,z,mh
    real Delta,chi,rvir,b,s,r,x,s0
    real rho,rmax

    rmax = 4

    integrate_rhotilde = 0

    Delta = bbps_delta(bbps_model)

    chi = rofz(z)
    rvir = (3*mh/4/pi/Delta/rhocrit(z))**(1./3.)     

    b = chi * sin(theta) / rvir 

    if(b>rmax) return

    s0 = sqrt(rmax**2-b**2) ! this is for a spherical cut at rmax
!    s0 = rmax               ! this is for a cylindrical cut at rmax

    n = int(s0/ds) + 1

    rho = 0
    do i=1,n
       s = (i-0.5) * ds
       x = sqrt(s**2+b**2)
       rho = rho + bbps_rhotilde(mh,z,x)
    enddo

    rho = rho * ds

    ! factor of two is because we only integrate half of it
    integrate_rhotilde = 2 * rho ! this thing has units Msun / kpc^3

    return 

  end function integrate_rhotilde

  real function integrate_ptilde_ell(ell,z,mh)

    implicit none

    integer i,n
    real ell,z,mh
    real Delta,chi,rvir,x
    real l2yl2

    Delta = bbps_delta(bbps_model)
    chi = rofz(z)
    rvir = (3*mh/4/pi/Delta/rhocrit(z))**(1./3.)     

    l2yl2 = 0
    x = 0
    do while(x<xmax)

       x     = x + dx/2
       l2yl2 = l2yl2 + x * bbps_ptilde(mh,z,x) * sin(ell*x*rvir/chi)
       x     = x + dx/2

    enddo
    l2yl2 = l2yl2**2

    integrate_ptilde_ell = l2yl2

    return 

  end function integrate_ptilde_ell

end module integrate_profiles
