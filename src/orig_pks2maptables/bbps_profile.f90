module bbps_profile

  use cosmology

  implicit none

  integer bbps_model

  real, dimension(3), parameter :: P0_Am     = (/ 18.100,  7.490, 20.7   /)
  real, dimension(3), parameter :: P0_alpham = (/  0.154,  0.226, -0.074 /)
  real, dimension(3), parameter :: P0_alphaz = (/ -0.758, -0.957, -0.743 /)
  
  real, dimension(3), parameter :: xc_Am     = (/  0.497,    0.710,  0.438 /)
  real, dimension(3), parameter :: xc_alpham = (/ -0.00865, -0.0833, 0.011 /)
  real, dimension(3), parameter :: xc_alphaz = (/  0.731,    0.853,  1.01  /)

  real, dimension(3), parameter :: beta_Am     = (/ 4.35,   4.19,  3.82   /)
  real, dimension(3), parameter :: beta_alpham = (/ 0.0393, 0.048, 0.0375 /)
  real, dimension(3), parameter :: beta_alphaz = (/ 0.415,  0.615, 0.535  /)

  real, dimension(3), parameter :: bbps_delta  = (/  200,  500,  200 /)
  real, dimension(3), parameter :: bbps_m2c    = (/ 0.26, 0.15, 0.26 /)

  real, parameter :: alpha = 1
  real, parameter :: gamma = -0.3

  real, dimension(3), parameter :: P0_Am_rho     = (/ 1.9e-5,  9.1e-05, 1.5e-4 /)
  real, dimension(3), parameter :: P0_alpham_rho = (/  0.29,  0.10, 0.14 /)
  real, dimension(3), parameter :: P0_alphaz_rho = (/ -0.67, -0.92, -1.32 /)

  real, dimension(3), parameter :: alpha_Am_rho     = (/  0.88,    0.69,  0.68 /)
  real, dimension(3), parameter :: alpha_alpham_rho = (/ -0.03, -0.02, -0.02 /)
  real, dimension(3), parameter :: alpha_alphaz_rho = (/  0.20,    0.27, 0.29  /)

  real, dimension(3), parameter :: beta_Am_rho     = (/ 3.83,   4.43,  6.4   /)
  real, dimension(3), parameter :: beta_alpham_rho = (/ 0.044, 0.0035, 0.028 /)
  real, dimension(3), parameter :: beta_alphaz_rho = (/ -0.025,  0.037, -0.25  /)

  real, parameter :: xc_rho = 0.5
  real, parameter :: gamma_rho = 0.2

  real, dimension(3), parameter :: P0_Am_dm     = (/ 2.4e-5,  2.4e-05, 1.5e-4 /)
  real, dimension(3), parameter :: P0_alpham_dm = (/  0.32,  0.30, 0.14 /)
  real, dimension(3), parameter :: P0_alphaz_dm = (/ 0.56, -0.84, -1.32 /)

  real, dimension(3), parameter :: alpha_Am_dm     = (/  0.85,    0.99,  0.68 /)
  real, dimension(3), parameter :: alpha_alpham_dm = (/ -0.03, -0.05, -0.02 /)
  real, dimension(3), parameter :: alpha_alphaz_dm = (/  0.18,    0.55, 0.29  /)

  real, dimension(3), parameter :: beta_Am_dm     = (/ 4.42,   3.40,  6.4   /)
  real, dimension(3), parameter :: beta_alpham_dm = (/ 0.072, 0.072, 0.028 /)
  real, dimension(3), parameter :: beta_alphaz_dm = (/ -0.22,  -0.81, -0.25  /)

  real, parameter :: xc_dm = 0.5
  real, parameter :: gamma_dm = 1.0


contains

  real function bbps_ptilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! The Compton y-parameter for a given halo is given by
    !   y = y0 * Mhalo * [Omegam(1+z)^3+Omegal] * int Ptilde * dltilde
    ! where
    !   Mhalo  = halo mass in Msun
    !   Ptilde = dimensionless pressure <--> Pressure = Ptilde*P_Delta
    !   ltilde = distance along los in units of virial radius rvir
    !   y0     = 3*sigma_T*(100 km/sec)^2/(8pi*me*c^2) in 1/Msun
    ! see 
    !   Battaglia et al. 2012, ApJ, 758, 75 
    ! section 4.1 for the precise definitions of P_Delta and Ptilde
    ! ----------------------------------------------------------------

    integer i

    real mh,z,x
    real A0,alpham,alphaz
    real P0,xc,beta

    i = bbps_model

    ! P0

    A0     = P0_Am(i)
    alpham = P0_alpham(i)
    alphaz = P0_alphaz(i)

    P0     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! xc

    A0     = xc_Am(i)
    alpham = xc_alpham(i)
    alphaz = xc_alphaz(i)

    xc     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! beta

    A0     = beta_Am(i)
    alpham = beta_alpham(i)
    alphaz = beta_alphaz(i)

    beta   = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! Pressure

    bbps_ptilde  = P0*(x/xc)**gamma*(1+(x/xc)**alpha)**-beta

    return

  end function bbps_ptilde

  real function bbps_rhotilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! The Compton y-parameter for a given halo is given y
    !   y = y0 * Mhalo * [Omegam(1+z)^3+Omegal] * int Ptilde * dltilde
    ! where
    !   Mhalo  = halo mass in Msun 
    !   Ptilde = dimensionless pressure <--> Pressure = Ptilde*P_Delta
    !   ltilde = distance along los in units of virial radius rvir
    !   y0     = 3*sigma_T*(100 km/sec)^2/(8pi*me*c^2) in 1/Msun
    ! see
    !   Battaglia in prep  2012, ApJ, 758, 75
    ! section 4.1 for the precise definitions of P_Delta and Ptilde
    ! ---------------------------------------------------------------- 

    integer i

    real mh,z,x
    real A0,alpham,alphaz
    real P0,beta,alpha_rho

    ! debug
    real r200c

    i = bbps_model

    ! P0

    A0     = P0_Am_rho(i)
    alpham = P0_alpham_rho(i)
    alphaz = P0_alphaz_rho(i)

     P0     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! xc

    A0     = alpha_Am_rho(i)
    alpham = alpha_alpham_rho(i)
    alphaz = alpha_alphaz_rho(i)

     alpha_rho     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! beta

    A0     = beta_Am_rho(i)
    alpham = beta_alpham_rho(i)
    alphaz = beta_alphaz_rho(i)

     beta   = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! gas density
    !bbps_rhotilde  = P0*(x/xc_rho)**(-gamma_rho)*(1+(x/xc_rho)**alpha_rho)**(-beta)
    bbps_rhotilde  = P0*(x/xc_rho)**(-1.0*gamma_rho)*(1+(x/xc_rho)**alpha_rho)**(-1.0*(beta-gamma_rho)/alpha_rho)
    ! SIS instead
!    r200c = (3.*mh/4./3.14159/200./rhocrit(z))**(1./3.) ! comoving Mpc
!    r200c = 1e3 * r200c / (1+z) ! proper kpc
!    bbps_rhotilde = 0.172 * mh / (4.*3.14159/3.*r200c**3) ! msun/kpc^3
!    bbps_rhotilde = bbps_rhotilde / 1e10  ! 1e10 h^2 msun/kpc^3

   return

  end function bbps_rhotilde

  real function bbps_dmtilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! The Compton y-parameter for a given halo is given y
    !   y = y0 * Mhalo * [Omegam(1+z)^3+Omegal] * int Ptilde * dltilde
    ! where
    !   Mhalo  = halo mass in Msun 
    !   Ptilde = dimensionless pressure <--> Pressure = Ptilde*P_Delta
    !   ltilde = distance along los in units of virial radius rvir
    !   y0     = 3*sigma_T*(100 km/sec)^2/(8pi*me*c^2) in 1/Msun
    ! see
    !   Battaglia in prep  2012, ApJ, 758, 75
    ! section 4.1 for the precise definitions of P_Delta and Ptilde
    ! ---------------------------------------------------------------- 

    integer i

    real mh,z,x
    real A0,alpham,alphaz
    real P0,beta,alpha_dm

    i = bbps_model

    ! P0

    A0     = P0_Am_dm(i)
    alpham = P0_alpham_dm(i)
    alphaz = P0_alphaz_dm(i)

    P0     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! xc

    A0     = alpha_Am_dm(i)
    alpham = alpha_alpham_dm(i)
    alphaz = alpha_alphaz_dm(i)

    alpha_dm     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! beta

    A0     = beta_Am_dm(i)
    alpham = beta_alpham_dm(i)
    alphaz = beta_alphaz_dm(i)

    beta   = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! Pressure

    bbps_dmtilde  = P0*(x/xc_dm)**-gamma_dm*(1+(x/xc_dm)**alpha_dm)**-beta

    return

  end function bbps_dmtilde


end module bbps_profile
