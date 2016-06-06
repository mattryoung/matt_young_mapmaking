module exclusion

  !-----------------------------------------------------------------------
  ! Routines for Lagrangian exclusion
  !-----------------------------------------------------------------------
  use mpivars
  integer, allocatable :: hoc(:,:,:)

contains 

  subroutine lagrangian_exclusion(x,y,z,vx,vy,vz,r,nin,nout)

!    implicit none

    integer nin, nout, nc, i, m, ix,iy,iz,imx,jmx,kmx,ll(nin),id(nin)
    real    x(nin),y(nin),z(nin),vx(nin),vy(nin),vz(nin),r(nin),&
            buffer(nin),f(nin)
    logical survived(nin), merged
    real    x1,x2,rmax,xmax

    parameter (nbr_mx=100000)
    real xlcl(nbr_mx),ylcl(nbr_mx),zlcl(nbr_mx),rlcl(nbr_mx)
    integer idlcl(nbr_mx)

    nc = 256

    allocate(hoc(0:nc-1,0:nc-1,0:nc-1))

    do m=1,nin
       id(m) = m
    enddo

    buffer = r
    call sort(nin,-r,id)
    r = buffer

    survived = .true.

    ! BUILD LINKED LISTS

    x1 = min(minval(x),minval(y),minval(z))
    x2 = max(maxval(x),maxval(y),maxval(z))

    ! CONVERT TO CELL SIZE LENGTH UNITS
    x = (x - x1)/(x2-x1) * nc
    y = (y - x1)/(x2-x1) * nc
    z = (z - x1)/(x2-x1) * nc
    r = (r     )/(x2-x1) * nc 

    hoc=0
    rmax = -1e10
    do mi=1,nin
       idi = id(mi)
       ix=int(x(idi))
       iy=int(y(idi))
       iz=int(z(idi))
       ix=min0(ix,nc-1)
       iy=min0(iy,nc-1)
       iz=min0(iz,nc-1)
       ll(mi)=hoc(ix,iy,iz)
       hoc(ix,iy,iz)=mi
    enddo

    ! LOOP OVER PEAKS, BEGINNING WITH MOST MASSIVE
    do mi=1,nin
       
       if(myid==0 .and. mod(mi,10000)==0) write(*,*) mi,nin,(100.*mi)/nin,'%'

       idi = id(mi)
       fi  = 0.0

       ! DON'T DO ANYTHING IF THIS PEAK HAS ALREADY BEEN ANNIHILATED
!       if(.not.survived(idi)) cycle

       xi  = x(idi)
       yi  = y(idi)
       zi  = z(idi)
       ri  = r(idi)

       i1 = max(0,   int(xi)-2*int(ri))
       i2 = min(nc-1,int(xi)+2*int(ri))

       j1 = max(0,   int(yi)-2*int(ri))
       j2 = min(nc-1,int(yi)+2*int(ri))

       k1 = max(0,   int(zi)-2*int(ri))
       k2 = min(nc-1,int(zi)+2*int(ri))

       nlcl = 0

       xlcl  = 0
       ylcl  = 0
       zlcl  = 0
       rlcl  = 0
       idlcl = 0

       ! FIRST MAKE LIST OF NEIGHBORS
       do i=i1,i2
          do j=j1,j2
             do k=k1,k2

                mj  = hoc(i,j,k)
                do while(mj>0)
                   idj = id(mj)
                   xj  =  x(idj)
                   yj  =  y(idj)
                   zj  =  z(idj)
                   rj  =  r(idj)
                   d = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
                   if (d < (ri + rj) .and. survived(idj) .and. ri >= rj .and. &
                        idj /= idi) then ! peaks i and j overlap
                      nlcl = nlcl + 1
                      if(nlcl > nbr_mx) stop
                      xlcl(nlcl)  = xj
                      ylcl(nlcl)  = yj
                      zlcl(nlcl)  = zj
                      rlcl(nlcl)  = rj
                      idlcl(nlcl) = idj
                   endif
                   mj = ll(mj)
                enddo

             enddo
          enddo
       enddo

       ! SORT NEIGHBORS BY MASS

       call sort(nlcl,-rlcl,idlcl)
       
       ! NOW GO DOWN THE LIST AND DO PAIRWISE EXCLUSION
       do mj=1,nlcl

          idj = idlcl(mj)

          if(survived(idj)) then
             xj = x(idj)
             yj = y(idj)
             zj = z(idj)
             rj = r(idj)
             
             call pairwise_exclusion(xi,yi,zi,ri,xj,yj,zj,rj,fi,merged)

             if(merged) survived(idj) = .false.

          endif

       enddo

       ! ASSIGN NEW VALUES TO PEAK i (displacements not re-calculated yet)

       x(idi) = xi
       y(idi) = yi
       z(idi) = zi
       r(idi) = r(idi) * (1-fi)**(1./3.)

    enddo

    ! FINAL CLEANUP
    nout = 0
    do i = 1,nin
       if(survived(i)) then
          nout = nout + 1
          x(nout) = x(i)
          y(nout) = y(i)
          z(nout) = z(i)
          r(nout) = r(i)
       endif
    enddo

    write(*,*) nin,nout

    ! CHANGE BACK TO PHYSICAL COORDINATES

    x = x1 + (x2-x1) / nc * x 
    y = x1 + (x2-x1) / nc * y
    z = x1 + (x2-x1) / nc * z
    r =      (x2-x1) / nc * r
       
    return

  end subroutine lagrangian_exclusion

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sort(nin,x,id)

    real x(nin),fid(nin)
    integer id(nin)

    fid = real(id)

    call sort2(nin,x,fid)

    id  = int(fid)

    return

  end subroutine sort

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pairwise_exclusion(xi,yi,zi,ri,xj,yj,zj,rj,f,merged)

    logical merged

    pi = 3.14159
    ftp = 4./3.*pi

    dij = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

    eta = 1.0
    merged = .false.
!    if(dij < ri ) then   ! half exclusion
    if(dij < ri+rj ) then ! full exclusion
       
       merged = .true.
!       if(dij**2 > ri**2 - rj**2) f = max(rj**3 / ri**3 * eta , f) 

    else
       write(*,*) ri,rj,dij

    endif

    return

  end subroutine pairwise_exclusion

  subroutine sphere_overlap(d,r1,r2,v1,v2)
    pi = 3.14159

    h1 = (r2 - r1 + d)*(r2 + r1 - d)/2./d
    h2 = (r1 - r2 + d)*(r1 + r2 - d)/2./d

    v1 = 1./3.*pi*h1**2*(3*r1-h1)
    v2 = 1./3.*pi*h2**2*(3*r2-h2)

    return 

  end subroutine sphere_overlap
         

end module exclusion
