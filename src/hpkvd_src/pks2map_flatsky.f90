program pks2map

  ! This code reads in a peak patch catalog and outputs a map.
  ! The contribution from each halo to a given pixel is stored 
  ! in a 3D table as a function of halo mass, redshift, and l.o.s.
  ! angle

  use cosmology
  use maptable
  use pksc
  use mpivars

  use pix_tools
  use fitstools
  use healpix_types
  use head_fits

  implicit none

  ! Filenames
  character *128 filein, fileout, tablefile
  character *4 outcode

  ! Compton y-parameter map
  real, allocatable :: field(:), fieldl(:)

  ! Map size and resolutions
  integer npix
  real fov,dpix,fov_cea,dp_cea

  ! Unit vectors, angles, and distances
  real xm,ym,zm,xp,yp,zp,xh,yh,zh,dx,dy,dz,x,y,z,d
  real theta,phi,theta0,phi0,theta1
  real dot,di,chi,redshift

  ! For FITS output
  logical simple,extend
  integer bitpix,naxis,naxes(2),group,ipix,unit,status

  ! External functions
  integer iargc, indlnb, i4arg
  real r4arg

  ! Map moments
  double precision mom1,mom2,mom3
  double precision hmom1,hmom2,hmom3,hmom1l,hmom2l,hmom3l

  ! Other variables
  integer i,j,k,jmin,jmax,kmin,kmax
  real theta_halo,phi_halo,dtheta_halo
  real thetamin,thetamax,phimin,phimax
  real thetap_cea,thetamin_cea,thetamax_cea
  real thetap,phip,dtheta,dphi,angle
  real y_halo,rh,ycur,ysigma
  real maxtheta,maxthetamass,maxthetachi,cut
  real mmin,chiview,zview,zcut
  real, parameter :: rmax = 4
  real mypi

  integer model

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)
  
  ! Usage check
  if(iargc()<4) then
     write(*,11) 
11   format('Usage: pks2map <filein> <fileout> <tablefile> <cut> [<npix> <fov> <chiview> <model>]')
     stop
  endif

  ! Get commandline
  call getarg(1,filein)
  call getarg(2,fileout)
  call getarg(3,tablefile)
  cut = r4arg(4,0.0)
  npix = i4arg(5,4096)
  fov = r4arg(6,2.0)
  chiview = r4arg(7,0)
  model = i4arg(8,1)

  mmin = 1e13

  testcase = .false.

  mypi = 4.*atan(1.0)

  ! Set field of view and resolution

  fov = fov / 360 * 2 * mypi ! radians
  fov_cea = 2 * sin(fov / 2)
  
  dpix = fov / npix
  dp_cea = fov_cea / npix

  if(myid==0) write(*,12) npix,npix,dpix,fov/2/mypi*360.
12 format(/,'Resolution is ',i0,'x',i0,' pixels (',1pe8.2,' radians)',/,&
            'Field of view is ',1pe8.2,' degrees')

  ! Allocate maps
  allocate(field(npix*npix))
  allocate(fieldl(npix*npix))

  ! Set background cosmology 
  omegam = 0.25
  omegab = 0.043
  omegal = 0.75
  h = 0.7
  sigma8 = 0.8
  ns = 0.96
  w = -1
  fb = omegab/omegam

  rho_0     = 2.775e11*omegam*0.7**2
  rhocrit_0 = 2.775e11*0.7**2  

  ! Create distance-redshift tables
  call set_rofztable
  call set_zofrtable

  zview = zofr(chiview)
  zcut = zofr(cut)
  if(myid==0) write(*,101) zview,zcut
101 format('zview ',f4.2,1x,f4.2)

  ! Load halos

  if(myid<10) then
     write(outcode,81) myid
  elseif(myid<100) then
     write(outcode,82) myid
  elseif(myid<1000) then
     write(outcode,83) myid
  else
     write(outcode,84) myid
  endif

  if(ntasks>1) filein=filein(1:indlnb(filein))//'.'//outcode

  call loadpksc(filein)

  mass = 4. / 3. * mypi * rho_0 * rth ** 3
  
  ! Load table

  call loadmaptable(tablefile)

  ! theta0 and phi0 are the center of fov
  ! here we set the center to lie along z-axis
  theta0 = 0
  phi0 = 0
  theta1 = mypi/2 - theta0
  
  ! Loop over halos

  fieldl=0
  maxtheta=0
  hmom1l=0
  hmom2l=0
  hmom3l=0
  do i=1,nhalo

     xh  = pos(3,i)
     yh  = pos(2,i) ! Note that x and z are switched
     zh  = pos(1,i)           

     chi = sqrt(xh**2+yh**2+zh**2)
     redshift=zofr(chi)

     mass(i) = sqrt(deltacrit(redshift)/bbps_delta(model)) * mass(i) ! Assuming SIS
     rh  = (3*mass(i)/4/mypi/delta_pksc/rhocrit(redshift))**(1./3.)

     xh = xh - chiview
     chi = sqrt(xh**2+yh**2+zh**2)

     if(chi<cut .or. mass(i)<mmin &
          .or. xh < 0) cycle

     dtheta_halo = rmax*asin(rh/chi)     

     ! Get unit vector

     x = xh / chi
     y = yh / chi
     z = zh / chi

     theta_halo = acos(z) 
     phi_halo   = asin(y/sin(theta_halo))
     theta_halo = theta_halo - mypi / 2

     thetamin = theta_halo - 1.1 * dtheta_halo
     thetamax = theta_halo + 1.1 * dtheta_halo

     phimin   = phi_halo - 1.1 * dtheta_halo
     phimax   = phi_halo + 1.1 * dtheta_halo

     if(thetamax < -fov/2 .or. thetamin > fov/2) cycle
     if(  phimax < -fov/2 .or.   phimin > fov/2) cycle

     if(thetamin <= -fov/2) thetamin = -fov/2 + 1e-10
     if(thetamax >=  fov/2) thetamax =  fov/2 - 1e-10

     if(  phimin <= -fov/2) phimin = -fov/2 + 1e-10
     if(  phimax >=  fov/2) phimax =  fov/2 - 1e-10
 
!     thetamin_cea = sin(thetamin)
!     thetamax_cea = sin(thetamax)

!     jmin = max(1,int((thetamin_cea + fov/2)/dp)+1)
!     jmax = min(int((thetamax_cea + fov/2)/dp)+1,npix)

     jmin = max(1,int((thetamin + fov/2)/dpix)+1)
     jmax = min(int((thetamax + fov/2)/dpix)+1,npix)

     kmin = max(1,int((phimin + fov/2)/dpix)+1)
     kmax = min(int((phimax + fov/2)/dpix)+1,npix)

     if(dtheta_halo > maxtheta) then
        maxtheta = dtheta_halo
        maxthetachi = chi
        maxthetamass = mass(i)
     endif
     
     do k=kmin,kmax
        phip   = -fov/2 + (k-0.5)*dpix
        do j=jmin,jmax           
!           thetap_cea = -fov_cea/2 + (j-0.5)*dp_cea
!           thetap = asin(thetap_cea)
           thetap   = -fov/2 + (j-0.5)*dpix

           dphi   = phip   - phi_halo
           dtheta = thetap - theta_halo

           angle = sqrt(dphi**2+dtheta**2)

           if(angle<dtheta_halo) then
              ycur =  interpolate_table(angle,redshift,mass(i))              
              fieldl(j+npix*(k-1)) = fieldl(j+npix*(k-1)) + ycur
              hmom1l = hmom1l + ycur
              hmom2l = hmom2l + ycur**2
              hmom3l = hmom3l + ycur**3
           endif

        enddo
     enddo
  enddo

  call mpi_allreduce(fieldl,field,npix**2,mpi_real,mpi_sum,mpi_comm_world,ierr)

  call mpi_allreduce(hmom1l,hmom1,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(hmom2l,hmom2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(hmom3l,hmom3,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  
  if(myid==0) then

  write(*,13) nhalo
13 format('Number of halos is ',i0)

  mom1 = sum(field**1) /npix / npix
  mom2 = sum(field**2) /npix / npix
  mom3 = sum(field**3) /npix / npix
  
  hmom1 = hmom1 / npix / npix
  hmom2 = hmom2 / npix / npix
  hmom3 = hmom3 / npix / npix

  write(*,14) zcut,mom1,hmom1,mom2,hmom2,mom3,hmom3,cut
14 format(8(1pe13.6,1x),'Moments')

  write(*,15) sqrt(sum((field-mom1)**2)/npix**2),mom1
15 format('RMS y-parameter fluctuation is ',1pe13.6,/,'Mean y-parameter is ',1pe8.2)
     
  write(*,16) maxtheta/2/mypi*360,maxthetamass,maxthetachi
16 format('Largest apparent halo is ',1pe8.2,' degrees, with mass ',1pe8.2,/,&
          ' at a distance of ',1pe8.2,'Mpc',/)

  if(testcase) field = field * 2.726e6 * 1.95 ! convert to |Delta t| in muK at 30 GHz

  open(unit=1,file=fileout(1:indlnb(fileout))//'.map',form='binary')
  write(1) npix,npix,fov,fov
  write(1) field
  close(1)

  if(ntasks==1) then

  ! Now output to FITS file
  
  simple   = .true.
  extend   = .true.
  bitpix   = 32
  naxis    = 2
  naxes(1) = npix
  naxes(2) = npix
  group    = 1
  ipix     = 1
  unit     = 11
  status   = 0

  ! Open file
  call ftinit(unit,fileout(1:indlnb(fileout))//'.fits',1,status)

  ! Write header
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! Write data
  call ftppre(unit,group,ipix,npix*npix,field,status)
  call ftpkye(unit,'SIZEX',fov,-1,'in radians',status)
  call ftpkye(unit,'SIZEY',fov,-1,'in radians',status)

  ! Close and free file
  call ftclos(unit,status)
  call ftfiou(unit,status)

  endif

  endif

  call mpi_finalize(ierr)

  stop

81 format('000',i1)
82 format( '00',i2)
83 format(  '0',i3)
84 format(      i4)

end program pks2map
