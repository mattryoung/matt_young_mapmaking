MODULE cosmology

  ! A MODULE THAT CONTAINS COSMOLOGICAL PARAMETERS AND FUNCTIONS

  implicit none
  
  real*8 tau_es,YHe,hubble0,ne0,thompson,clight,dtau0,dtau
  integer NHe

  real omegam,omegal,omegab,h,w,ns,sigma8
  real obhfactor,omfactor
  real taufactor
  real tcmb0,zinit,tcti,c,h0,rinit
  
  parameter(omegam=0.27,omegal=1.-omegam,omegab=0.042,h=0.72)
  parameter(w=-1,ns=1,sigma8=0.8)

  parameter(obhfactor=omegab*h/0.03,omfactor=1./(omegam/0.27)**0.5)
  parameter(taufactor=3.3e-3*obhfactor*omfactor)

  parameter(tcmb0=2.726)
  parameter(zinit=100000.,rinit=14.23)
  parameter(h0=1e5*h,c=3.e5)  ! h0 in [km/s/Gpc] and c in [km/s] !!
  parameter(tcti=2.*c/h0/omegam**0.5/(1.+zinit)**1.5/3.086e27) ! tcti = 3cti in Gpc

  parameter(clight=3e10,thompson=6.65e-25,YHe=0.25)
  parameter(hubble0=100./3.086e19*h,ne0=omegab*h**2*1.88e-29/1.67e-24)
  parameter(dtau0=thompson*ne0*clight/hubble0)

  real zofrtable_min,zofrtable_max,dzofrtable
  integer nzofrtable  

  real, allocatable:: zofrtable(:)

  real mmin ! MINIMUM HALO MASS
  
CONTAINS
  
  !=======================================================================

  REAL FUNCTION dtaudz(x,z)
    
    implicit none
    
    real x,z
    
    dtaudz=x*sqrt(1.+z)
    
  END FUNCTION dtaudz
  
  !=======================================================================

  subroutine set_zofrtable
    
    implicit none

    real z
    integer i,j

    zofrtable_min = 0.
    zofrtable_max = 16.
    nzofrtable = 100000
    dzofrtable = (zofrtable_max-zofrtable_min)/(nzofrtable-1)

    allocate(zofrtable(nzofrtable))

    zofrtable(1)=0.

    do i=2,nzofrtable
       z=zofrtable(i-1)
       if(z>1.e6) then
          do j=i,nzofrtable
             zofrtable(j)=z
          enddo
          return
       endif
       zofrtable(i)=zofrtable(i-1)+dzdr(z)*dzofrtable
    enddo

    return

  end subroutine set_zofrtable

  !=======================================================================

  real function dzdr(z)

    implicit none

    real z

    dzdr = h0/c*sqrt(omegam*(1.+z)**3+1.-omegam)

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

  real FUNCTION zofr_old(r)
    
    real r
    real delta
    
    !     dr = c * dz/H(z) = c / H0 / sqrt(omegam) * dz/(1+z)^1.5
    !     r - ri = c / H0 / sqrt(omegam) * int[z..zi] dz/(1+z)^1.5
    !     r - ri = c * 2 / H0 / sqrt(omegam) * [(1+z)^-0.5 - (1+zi)^-0.5]
    !     z = [2 * (r-ri) * H0 / c * sqrt(omegam) + (1+zi)^-0.5]^-2 - 1
    
    if(r.ge.rinit) then
       zofr_old = zinit
       return
    endif
    
    delta = (rinit-r)*3.086e27*h0/c/2.*sqrt(omegam)+(1.+zinit)**(-0.5)
    
    zofr_old = delta**(-2) - 1.
    
    return
    
  END FUNCTION zofr_old
  
END MODULE cosmology
