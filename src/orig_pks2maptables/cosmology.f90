module cosmology

  implicit none

  real omegam,omegab,omegal,h,ns,sigma8,w,h0,fb ! h0 is in km/s/Mpc
  real rho_0, rhocrit_0, y0, tau0, mue

  double precision pi
  parameter (pi = 3.14159265359)

  double precision sigmat,c,me,mp,msun, Mpc, kpc ! all in cgs
  double precision Mpc2km, H100, ckms            ! H100 is 100 km/s/Mpc in 1/s
  parameter (sigmat = 6.65246d-25  ,&
                  c = 2.99792458d10,&
                 me = 9.109383d-28 ,&
                mue = 1.136        ,&
                 mp = 1.67e-24     ,&
               msun = 1.9891d33    ,&
             Mpc2km = 3.08567758d19,&
                Mpc = 3.08567758d24,&
                kpc = 3.08567758d21,&
               H100 = 1e2/Mpc2km   ,&
               ckms = c/1e5)

  double precision ypermass ! this is the quantity
                            ! 3*sigmat*ke*(100 km/sec/Mpc)^2/(8*pi*me*c^2)
                            ! in units of 1/Msun
  parameter (ypermass = 3*sigmat*H100**2/8./pi/me/c**2*msun/1.932)

  double precision taupersigma ! this is the quantity 
                               ! sigmat/mu/mp
                               ! in units of 1/(1e10 Msun/kpc^3)/Mpc
  parameter (taupersigma = sigmat/mue/mp*msun/kpc**3*Mpc*1e10)

  double precision ylpermsun ! this is the quantity
                             ! 3*sigmat*ke*(100 km/sec/Mpc)^2/4/(me*c^2)
                             ! in units of 1/Msun
  parameter (ylpermsun = 3*sigmat*H100**2/4./me/c**2*msun/1.932)

  integer nzofrtable,nrofztable
  real zofrtable_min,zofrtable_max,dzofrtable
  real rofztable_min,rofztable_max,drofztable
  real, allocatable :: rofztable(:),zofrtable(:)

  logical testcase 

contains

  !=======================================================================

  subroutine set_rofztable
    
    implicit none
    
    real r,z,dz
    integer i,j
    
    rofztable_min = 0
    rofztable_max = 6.
    nrofztable = 10000
    drofztable = (rofztable_max-rofztable_min)/(nrofztable-1)

    dz = drofztable

    allocate(rofztable(nrofztable))

    rofztable(1)=0.

    do i=2,nrofztable
       z=(i-1)*dz+rofztable_min
       rofztable(i)=rofztable(i-1)+drdz(z)*dz
    enddo

    return

  end subroutine set_rofztable

  !=======================================================================  

  subroutine set_zofrtable
    
    implicit none
    
    real z,dr
    integer i,j
    
    zofrtable_min = 0.
    zofrtable_max = 1.3e4
    nzofrtable = 10000
    dzofrtable = (zofrtable_max-zofrtable_min)/(nzofrtable-1)

    dr=dzofrtable

    allocate(zofrtable(nzofrtable))

    zofrtable(1)=0.

    do i=2,nzofrtable
       z=zofrtable(i-1)
       zofrtable(i)=zofrtable(i-1)+dzdr(z)*dr
    enddo

    return

  end subroutine set_zofrtable

  !=======================================================================

  real function drdz(z)

    implicit none

    real z

    drdz = 1 / dzdr(z)

    return

  end function drdz

  !=======================================================================

  real function dzdr(z)

    implicit none

    real z

    dzdr = h*100/ckms*sqrt(omegam*(1.+z)**3+1.-omegam) ! dzdr is in 1/Mpc

    return
    
  end function dzdr

  !=======================================================================

  real function zofr(r)

    real r
    real fint
    integer itab

    itab = ((r-zofrtable_min)/dzofrtable)+1

    if(itab<1) then
       zofr = zofrtable(1)
       return
    elseif(itab>=nzofrtable) then
       zofr = zofrtable(nzofrtable)
       return
    endif

    fint=(r-real(itab-1)*dzofrtable-zofrtable_min)/dzofrtable

    zofr=(1-fint)*zofrtable(itab)+fint*zofrtable(itab+1)

  end function zofr

  !=======================================================================

  real function rofz(z)

    real z
    real fint
    integer itab

    itab = ((z-rofztable_min)/drofztable)+1

    if(itab<1) then
       rofz = rofztable(1)
       return
    elseif(itab>=nrofztable) then
       rofz = rofztable(nrofztable)
       return
    endif

    fint=(z-real(itab-1)*drofztable-rofztable_min)/drofztable

    rofz=(1-fint)*rofztable(itab)+fint*rofztable(itab+1)

  end function rofz

  !=======================================================================

  real function rhocrit(z)

    ! This is the *comoving* critical density as a function of redshift

    real z

    rhocrit = rhocrit_0 * (omegam*(1+z)**3+omegal) / (1+z)**3

    return

  end function rhocrit

  !=======================================================================

  real function deltacrit(z)

    ! This is the post-spherical-collapse overdensity for an isolated 
    ! uniform sphere in virial equilibrium in units of the *critical*
    ! density; see Bryan & Norman (1997)

    real z,x,omegaz

    omegaz = omegam*(1+z)**3/(omegam*(1+z)**3+1-omegam)
    x = omegaz - 1

    deltacrit = 18*pi**2 + 82 * x - 39 * x**2

    return

  end function deltacrit

end module cosmology

