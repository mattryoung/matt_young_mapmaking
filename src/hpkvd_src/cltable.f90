module cltable

  use cosmology
  use bbps_profile
  use integrate_profiles

  implicit none

  ! Table bounds and interval
  real     zt,   zmint,   zmaxt
  real    mht,  mhmint,  mhmaxt
  real   ellt, ellmint, ellmaxt

  ! Table sizes and intervals
  integer nellt, nzt, nmht  
  real    dellt, dzt, dmht

  ! Table
  real, allocatable :: table(:,:,:)

contains

  subroutine makecltable(fileout,model)

    implicit none 

    character *128 fileout
    integer model

    ! lyl(M,z)
    real elly,l2yl2
    real yl0

    ! Other variables
    integer i,j,k
    real z,theta,Delta,zfac

    ! Create distance-redshift tables
    call set_rofztable
    call set_zofrtable
  
    bbps_model = model
    Delta = bbps_delta(bbps_model)
    
    yl0 = fb * Delta * h**2 * ylpermsun
 
    ! Setup table bounds
    
    ! Redshift
    nzt   = 100
    zmint = 1e2
    zmaxt = 7e3
    
    dzt = (log(zmaxt)-log(zmint))/(nzt-1)
    
    ! Mass
    nmht   = 100
    mhmint = 1e13
    mhmaxt = 1e16
    
    dmht = (log(mhmaxt)-log(mhmint))/(nmht-1)

    ! multipole
    nellt   = 200
    ellmint = 10
    ellmaxt = 1e5
    
    dellt = (log(ellmaxt)-log(ellmint))/(nellt-1)

    allocate(table(nellt,nmht,nzt))
    
    ! Calculate tabulated values and output
    
    open(unit=1,file=fileout,form='binary')
    write(1) nzt,nmht,nellt,zmint,zmaxt,mhmint,mhmaxt,ellmint,ellmaxt
    
    do i=1,nzt
       zt = exp( log(zmint) + (i-1) * dzt )
       zfac = (1+zt)**3
       write(*,*) i,nzt

       do j=1,nmht
          mht = exp( log(mhmint) + (j-1) * dmht )

          do k=1,nellt
             ellt = exp( log(ellmint) + (k-1) * dellt )

             l2yl2 = integrate_ptilde_ell(ellt,zt,mht)
             l2yl2 = yl0 * zfac * mht * l2yl2
             
             table(k,j,i) = l2yl2 ! This is to preserve C-like array order in
                                  ! header without thrashing the cache

          enddo
       enddo
    enddo

    write(1) table
    
    close(1)
    
  end subroutine makecltable
  
  !=======================================================================

  subroutine loadcltable(filein)

    implicit none

    character *128 filein

    open(unit=1,file=filein,form='binary')

    read(1) nzt,nmht,nellt,zmint,zmaxt,mhmint,mhmaxt,ellmint,ellmaxt

    allocate(table(nellt,nmht,nzt))

    read(1) table

    close(1)

    dzt   = (log(  zmaxt) - log(  zmint)) / (  nzt - 1)
    dellt = (log(ellmaxt) - log(ellmint)) / (nellt - 1)
    dmht  = (log( mhmaxt) - log( mhmint)) / ( nmht - 1)

    return

  end subroutine loadcltable

  !=======================================================================

  real function interpolate_cltable(ell,z,mh)

    implicit none

    real z,ell,mh
    real fell,fz,fm
    integer iell,iz,imh

    interpolate_cltable = 0
    if(ell>ellmaxt .or. z>zmaxt) return

    if(ell < ellmint) ell = ellmint + 1e-5
    if(  z <   zmint)   z =   zmint + 1e-5
    if( mh <  mhmint)  mh =  mhmint + 1e-5
    if( mh >  mhmaxt)  mh =  mhmaxt - 1e-5

    iell = int( ( log(ell) - log(ellmint) ) / dellt ) + 1
    iz   = int( ( log(  z) - log(  zmint) ) / dzt   ) + 1
    imh  = int( ( log( mh) - log( mhmint) ) / dmht  ) + 1

    fell = log(ell) - ( log(ellmint) + (iell - 1) * dellt )
    fz   = log(  z) - ( log(  zmint) + (  iz - 1) *   dzt )
    fm   = log( mh) - ( log( mhmint) + ( imh - 1) *  dmht )
    
    interpolate_cltable = &
         table(iell  ,imh  ,iz  ) * (1-fell) * (1-fz) * (1-fm) + &
         table(iell+1,imh  ,iz  ) * (  fell) * (1-fz) * (1-fm) + &
         table(iell  ,imh+1,iz  ) * (1-fell) * (  fz) * (1-fm) + &
         table(iell+1,imh+1,iz  ) * (  fell) * (  fz) * (1-fm) + &
         table(iell  ,imh  ,iz+1) * (1-fell) * (1-fz) * (  fm) + &
         table(iell+1,imh  ,iz+1) * (  fell) * (1-fz) * (  fm) + &
         table(iell  ,imh+1,iz+1) * (1-fell) * (  fz) * (  fm) + &
         table(iell+1,imh+1,iz+1) * (  fell) * (  fz) * (  fm) 
    
    interpolate_cltable = interpolate_cltable 
            
  end function interpolate_cltable

end module cltable
