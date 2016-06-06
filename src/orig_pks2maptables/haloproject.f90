module haloproject

  use cosmology

  use mpivars
  use flatskyvars
  use healpixvars
  use maptable

  use bbps_profile

  use pix_tools
  use healpix_types
  use head_fits

  implicit none

  real, parameter :: rmax = 4
  double precision theta_halo,phi_halo,dtheta_halo
  real angle

  ! Halo angular properties
  real maxtheta,maxthetamass,maxthetachih  
  real maxthetal,maxthetamassl,maxthetachihl  

  ! Map moments
  double precision fsmom1,fsmom2,fsmom3
  double precision hpmom1,hpmom2,hpmom3
  real curval

  ! Weird stuff
  real localpi

contains

!=======================================================================

  subroutine projecthalo_healpix(x,y,z,r,chi,m,redshift,v,profile)    

    implicit none

    real x,y,z,r,chi,m,redshift,v
    integer i,j,profile

    localpi = 4*atan(1.0)

    dtheta_halo = rmax*asin(r/chi)
    
    ! Get unit vector
    
    hpvec1(1) = z / chi ! was x,y,z
    hpvec1(2) = y / chi
    hpvec1(3) = x / chi   

    call query_disc(hpnside,hpvec1,dtheta_halo,hplist,nlist)
    
!    write(*,*) "theta =0, z=0, m=1e15 halo value = ", interpolate_table(0.0,0.0,1e15,2) 
!    stop
    
    do i=0,nlist-1

       j=hplist(i)
       call pix2vec_ring(hpnside,j,hpvec2)
       call angdist(hpvec1,hpvec2,hpangle)
       
       angle = hpangle

       if(j>hpnpix-1.or.j<0) write(*,*) 'error',i,j,hpnpix,nlist
       
       if(angle<dtheta_halo) then
          curval =  v * interpolate_table(angle,redshift,m,profile)              
          hpmapl(j,1) = hpmapl(j,1) + curval
       endif

    enddo

    return

  end subroutine projecthalo_healpix

!=======================================================================

  subroutine projecthalo_flatsky(x,y,z,r,chi,m,redshift,v,profile)

    implicit none

    real x,y,z,r,chi,m,redshift,v

    integer i,j,k,jmin,jmax,kmin,kmax,profile

    localpi = 4*atan(1.0)

    dtheta_halo = rmax*asin(r/chi)
    
    ! Get unit vector
    
    x = x / chi
    y = y / chi
    z = z / chi

    if(z<0.0) return
!    theta_halo = acos(z) 
!    phi_halo   = asin(y/sin(theta_halo))
!    theta_halo = theta_halo - localpi / 2
    theta_halo = asin(x) 
    phi_halo   = asin(y)

    
    thetamin = theta_halo - 1.1 * dtheta_halo
    thetamax = theta_halo + 1.1 * dtheta_halo
    
    phimin   = phi_halo - 1.1 * dtheta_halo
    phimax   = phi_halo + 1.1 * dtheta_halo
    
    if(thetamax < -fov/2 .or. thetamin > fov/2) return
    if(  phimax < -fov/2 .or.   phimin > fov/2) return
    
    if(thetamin <= -fov/2) thetamin = -fov/2 + 1e-10
    if(thetamax >=  fov/2) thetamax =  fov/2 - 1e-10
    
    if(  phimin <= -fov/2) phimin = -fov/2 + 1e-10
    if(  phimax >=  fov/2) phimax =  fov/2 - 1e-10
    
    jmin = max(1,int((thetamin + fov/2)/dpix)+1)
    jmax = min(int((thetamax + fov/2)/dpix)+1,npix)
    
    kmin = max(1,int((phimin + fov/2)/dpix)+1)
    kmax = min(int((phimax + fov/2)/dpix)+1,npix)
    
    if(dtheta_halo > maxtheta) then
       maxtheta = dtheta_halo
       maxthetachih = chi
       maxthetamass = m
    endif
    
    do k=kmin,kmax
       phip   = -fov/2 + (k-0.5)*dpix
       do j=jmin,jmax           
          thetap   = -fov/2 + (j-0.5)*dpix
          
          dphi   = phip   - phi_halo
          dtheta = thetap - theta_halo
          
          angle = sqrt(dphi**2+dtheta**2)
          
          if(angle<dtheta_halo) then
             curval =  v*interpolate_table(angle,redshift,m,profile)              
             fsmapl(j+npix*(k-1)) = fsmapl(j+npix*(k-1)) + curval
          endif
          
       enddo
    enddo
    
    return 
    
  end subroutine projecthalo_flatsky
  
end module haloproject
