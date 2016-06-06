module pksc  

  use cosmology
  use mpivars 

  implicit none

  real delta_pksc
  parameter(delta_pksc=200)

  ! Halo arrays
  real, allocatable :: posxyz(:,:),velxyz(:,:),rth(:),mass(:),vrad(:)
  real rhi
  integer nhalo, nhalol 
  real RTHmaxtot
  integer offset_num_floats

contains

  !=======================================================================

  subroutine loadpksc(filename,profile)
    
    implicit none
    
    integer i,j,idum,m,profile

    character *128 filename

    if(testcase) then

       nhalo = 1
       
       allocate(posxyz(3,nhalo))
       allocate(rth(nhalo))
       allocate(mass(nhalo))
       
       posxyz(1,1) = 0
       posxyz(2,1) = 0
       posxyz(3,1) = rofz(0.05)
             
       return
       
    endif

    open(unit=1,file=filename,form='binary')
    read(1) nhalo, RTHmaxtot
    close(1)
    
    !Read in only the appropriate nhalo/ntasks, not whole catalogue
    nhalol = int((nhalo-1)/ntasks)+1          ! even number per processor
    nhalo  = min(nhalol, nhalo - nhalol*myid) ! last processor may not have 
                                              ! same number of halos
    allocate(posxyz(3,nhalo))
    allocate(velxyz(3,nhalo))
    allocate(rth(nhalo))
    allocate(mass(nhalo))
    allocate(vrad(nhalo))

    offset_num_floats = 7*nhalo*myid
    open(unit=1,file=filename,form='binary')
    read(1) idum, idum
    read(1) (idum,i=1,offset_num_floats) !read into dummy all halos previous
    read(1) ((posxyz(j,i),j=1,3),&
             (velxyz(j,i),j=1,3),&
             rth(i),i=1,nhalo)

    close(1)

    if (profile==2) then
       vrad = (posxyz(:,1)*velxyz(:,1) + posxyz(:,2)*velxyz(:,2)+ posxyz(:,3)*velxyz(:,3))
       vrad = vrad/(posxyz(:,1)**2 + posxyz(:,2)**2 + posxyz(:,3)**2)**(1./2)
       write(*,*) "VRAD min max = ", minval(vrad), maxval(vrad)
       vrad = -vrad / 3.e5 !negative as ksz is -v/c*tau
!       vrad = 1.0
    elseif (profile==1) then
       vrad = 1.0 
    endif
    
!    if(myid==0) then
!       write(*,*) "min max x,y,z,rth"
!       write(*,*) minval(posxyz(1,:)),maxval(posxyz(1,:))
!       write(*,*) minval(posxyz(2,:)),maxval(posxyz(2,:))
!       write(*,*) minval(posxyz(3,:)),maxval(posxyz(3,:))
!       write(*,*) minval(rth),maxval(rth)
!    endif

    return

  end subroutine loadpksc

end module pksc

