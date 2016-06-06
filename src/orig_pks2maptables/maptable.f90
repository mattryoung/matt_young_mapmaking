module maptable

  use cosmology
  use bbps_profile
  use integrate_profiles

  implicit none

  ! Table bounds and interval
  real   chit,   chimint,   chimaxt
  real    mht,    mhmint,    mhmaxt
  real     rt,     rmint,     rmaxt

  ! Table sizes and intervals
  integer nrt, nchit, nmht, nprofiles
  real    drt, dchit, dmht

  parameter (nprofiles=2)

  ! Table
  real, allocatable :: table(:,:,:,:)

contains

  subroutine makemaptable(fileout,model)

    implicit none 

    character *128 fileout
    integer model

    ! Compton y-parameter
    real y

    ! Thomson scattering optical depth
    real tau

    ! Other variables
    integer i,j,k
    real z,theta,Delta,zfacy,zfact,m2rfac,mht3
    real r200c

    ! Debugging
    real tau_mean, tau_mean_expected, thetap, dtheta

    ! Create distance-redshift tables
    call set_rofztable
    call set_zofrtable
  
    bbps_model = model
    Delta = bbps_delta(bbps_model)
    
    m2rfac = (3./4./pi/Delta)**(1./3.) 

    y0     = fb * Delta * h**2 * ypermass  ! in 1/Msun
    tau0   = taupersigma * h**2 * m2rfac   ! 1/(1e10 h^2*Msun/kpc^3)/Mpc
 
    ! Setup table bounds
    
    ! Distance
    nchit   = 100
    chimint = 1e2
    chimaxt = 7e3
    
    dchit = (log(chimaxt)-log(chimint))/(nchit-1)
    
    ! Mass
    nmht   = 100
    mhmint = 1e13
    mhmaxt = 1e16
    
    dmht = (log(mhmaxt)-log(mhmint))/(nmht-1)

    ! Radius
    nrt   = 200 !200
    rmint = 1e-2
    rmaxt = 10 !10
    
    drt = (log(rmaxt)-log(rmint))/(nrt-1)

    allocate(table(nprofiles,nrt,nmht,nchit))
    
    ! Calculate tabulated values and output
    
    open(unit=1,file=fileout,form='binary')
    write(1) nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt
    
    do i=1,nchit
       chit  = exp( log(chimint) + (i-1) * dchit )
       z  = zofr(chit)
       zfacy = (omegam*(1+z)**3+omegal)
       zfact = 1./(1+z)/rhocrit(z)**(1./3.) ! Mpc/Msun^(1/3)

       do j=1,nmht
          mht  = exp( log(mhmint) + (j-1) * dmht )
          mht3 = mht**(1./3.)

          tau_mean_expected = 1e-9 * (mht/1e15) / (chit/1e3)**2 
          tau_mean = 0
          thetap = 0.
          do k=1,nrt
             rt = exp( log(rmint) + (k-1) * drt )
             theta = atan(rt/chit)
             dtheta = theta-thetap
             y   =   y0 * mht  * zfacy * integrate_ptilde(  theta,z,mht)

             ! [tau0]  = 1/(1e10 Msun/kpc^3)/Mpc
             ! [zfact] = Mpc/Msun^(1/3)
             ! [mht3]  = Msun(1/3)
             ! [integrate_rhotilde] = 1e10 Msun/kpc^3
             tau = tau0 * mht3 * zfact * integrate_rhotilde(theta,z,mht)

             table(1,k,j,i) = y 
             table(2,k,j,i) = tau

             tau_mean = tau_mean + dtheta * theta * tau

             thetap = theta
          enddo
          tau_mean = tau_mean / 2.
       enddo
    enddo

    write(1) table
    
    close(1)
    
  end subroutine makemaptable
  
  !=======================================================================

  subroutine loadmaptable(filein)

    implicit none

    character *128 filein

    open(unit=1,file=filein,form='binary')

    read(1) nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt

    allocate(table(nprofiles,nrt,nmht,nchit))

    read(1) table

    close(1)

    dchit = (log(chimaxt) - log(chimint)) / (nchit - 1)
    drt   = (  log(rmaxt) -   log(rmint)) / (nrt - 1)
    dmht  = ( log(mhmaxt) -  log(mhmint)) / (nmht - 1)

    return

  end subroutine loadmaptable

  !=======================================================================

  real function interpolate_table(theta,z,mh,profile)

    implicit none

    real theta,z,chi,r,mh
    real fr,fc,fm
    integer ir,ichi,imh,profile

    chi = rofz(z)
    r   = chi * tan(theta)

    interpolate_table = 0
    if(r>rmaxt .or. chi>chimaxt) return

    if(   r < rmint   )   r =   rmint + 1e-5
    if( chi < chimint ) chi = chimint + 1e-5
    if(  mh < mhmint  )  mh =  mhmint + 1e-5

    if(  mh > mhmaxt  )  mh =  mhmaxt - 1e-5

    ir   = int( (   log(r) -   log(rmint) ) / drt   ) + 1
    ichi = int( ( log(chi) - log(chimint) ) / dchit ) + 1
    imh  = int( (  log(mh) -  log(mhmint) ) / dmht  ) + 1

    fr = log(r)   - (   log(rmint) + (  ir - 1) *   drt )
    fc = log(chi) - ( log(chimint) + (ichi - 1) * dchit )
    fm = log(mh)  - (  log(mhmint) + ( imh - 1) *  dmht )
    
    interpolate_table = &
         table(profile,ir  ,imh  ,ichi  ) * (1-fr) * (1-fc) * (1-fm) + &
         table(profile,ir+1,imh  ,ichi  ) * (  fr) * (1-fc) * (1-fm) + &
         table(profile,ir  ,imh+1,ichi  ) * (1-fr) * (  fc) * (1-fm) + &
         table(profile,ir+1,imh+1,ichi  ) * (  fr) * (  fc) * (1-fm) + &
         table(profile,ir  ,imh  ,ichi+1) * (1-fr) * (1-fc) * (  fm) + &
         table(profile,ir+1,imh  ,ichi+1) * (  fr) * (1-fc) * (  fm) + &
         table(profile,ir  ,imh+1,ichi+1) * (1-fr) * (  fc) * (  fm) + &
         table(profile,ir+1,imh+1,ichi+1) * (  fr) * (  fc) * (  fm) 
    
    interpolate_table = interpolate_table 
            
  end function interpolate_table

end module maptable
