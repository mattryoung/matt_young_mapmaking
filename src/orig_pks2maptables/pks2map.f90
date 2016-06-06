program pks2map

  ! This code reads in a peak patch catalog and outputs a map.
  ! The contribution from each halo to a given pixel is stored 
  ! in a 3D table as a function of halo mass, redshift, and l.o.s.
  ! angle

  use cosmology

  use pksc
  use haloproject

  use mpivars
  use fitsvars
  use healpixvars
  use flatskyvars

  use fitstools
  use healpix_types
  use head_fits

  use textlib
  
  use random

  implicit none

  ! Filenames
  character *128 filein, fileout, fileout_bn,fileout_fs,fileout_hp,&
       fileout_dp,tablefile
  character *4 outcode
  integer ilast

  ! Halo properties
  real xh,yh,zh,mh,rh,redshifth,chih,vh
  real xmmax,xmmaxl, ymmax,ymmaxl, zmmax,zmmaxl, rmmax

  ! dPdy
  integer, parameter :: ndpdy = 10000
  real,    parameter :: dpdy_max = 1e-4,dpdy_dy=dpdy_max/ndpdy
  real dpdy(ndpdy)

  ! Other variables
  integer i,j,k,jmin,jmax,kmin,kmax,m

  real mmin,chihview,zview,cut_low,cut_high,zmin,zmax
  real mypi
  real RTHmaxl,RTHmax,rmaxtheta,rmaxphi
  real muran, phiran
  integer nhalotot
  integer model, centre, scramble, profile

  ! For random number generation
  integer,dimension( IRandNumSize ) :: seed
  real(kind=8) :: randnum

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)
  
  ! Usage check
  if(iargc()<3) then
     if(myid==0) write(*,11) 
11   format('Usage: pks2map <filein> <fileout> <tablefile> [<zmin> <nside> <scramble> <center> <npix> <fov> <zmax> <chihview> <model> <profile>]')
     call mpi_finalize(ierr)
     stop
  endif

  ! Get commandline
  call      getarg(1,filein)
  call      getarg(2,fileout)
  call      getarg(3,tablefile)
  zmin     = r4arg(4,0.0)
  hpnside  = i4arg(5,4096)
  scramble = i4arg(6,0)
  centre   = i4arg(7,0) 
  npix     = i4arg(8,16)
  fov      = r4arg(9,10.0)
  zmax     = r4arg(10,1.25)
  chihview = r4arg(11,0.0)
  model    = i4arg(12,1)
  profile  = i4arg(13,1)

  mmin = 1e13

  testcase = .false.

  mypi = 4.*atan(1.0)

  ! Set field of view and resolution

  fov     = fov / 360 * 2 * mypi ! radians
  dpix   = fov     / npix

  if(myid==0) write(*,12) npix,npix,dpix,fov/2/mypi*360.
12 format(/,/,&
          'Resolution:         ',i0,' x ',i0,' pixels (',1pe8.2,' radians)',/,&
          'Field of view:      ',1pe8.2,' degrees')

  ! Healpix map resolution

  hpnpix = nside2npix(hpnside)

  if(myid==0) write(*,13) hpnside,hpnpix
13 format('Healpix Nside:      ',i0,/,&
          'Number of pixels:   ',i0)

  ! Allocate maps

  allocate(fsmap(npix*npix))
  allocate(fsmapl(npix*npix))

  allocate(hpmap(0:hpnpix-1,1:1))
  allocate(hpmapl(0:hpnpix-1,1:1))
  allocate(hplist(0:hpnpix-1))

  ! Set background cosmology 
  omegam = 0.25
  omegab = 0.043
  omegal = 0.75
  h = 0.7
  sigma8 = 0.8
  ns = 0.96
  w = -1
  fb = omegab/omegam

  rho_0     = 2.775e11*omegam*h**2
  rhocrit_0 = 2.775e11*h**2  

  ! Create distance-redshift tables
  call set_rofztable
  call set_zofrtable

  cut_low  = rofz(zmin)
  cut_high = rofz(zmax)

  if(myid==0) write(*,14) cut_low/1e3,cut_high/1e3
14 format('cut_low, cut_high:  ',f4.2,1x,f4.2,' Gpc')

  ! Load halos
  call loadpksc(filein,profile)

  ! Cut all the halos not within range
  m=0
  do i=1,nhalo 
     xh  = posxyz(3,i)
     yh  = posxyz(2,i) ! Note that x and z are switched
     zh  = posxyz(1,i)           
     xh = xh - chihview

     chih = sqrt(xh**2+yh**2+zh**2)
     if(chih>cut_high.or.chih<cut_low) cycle

     m = m+1
     posxyz(3,m) = posxyz(3,i)
     posxyz(2,m) = posxyz(2,i)
     posxyz(1,m) = posxyz(1,i)
     rth(m)      = rth(i)
  enddo
  nhalo = m

  if(scramble==1) then
     ! Scramble theta and phi
     seed = 13579 * (myid+1)
     do i=1,nhalo 
        xh   = posxyz(1,i)
        yh   = posxyz(2,i)
        zh   = posxyz(3,i)
        chih = sqrt(xh**2+yh**2+zh**2)

        call ranf(seed,randnum)
        muran  = 2 * randnum - 1
        call ranf(seed,randnum)
        phiran = 2 * 3.14159 * randnum

        rh = sqrt(1-muran**2) * chih

        xh = rh    * cos(phiran)
        yh = rh    * sin(phiran)
        zh = muran * chih

        posxyz(1,i) = xh
        posxyz(2,i) = yh
        posxyz(3,i) = zh

     enddo
  endif

  if(centre==0) goto 111 ! if don't want a centered map on most massive halo

  !find largest remaining halo across all processors
  RTHmaxl = maxval(rth(1:nhalo))
  call mpi_allreduce(RTHmaxl,RTHmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)

  ! find x,y,z position amd angles of largest halo
  xmmaxl = -1e5
  ymmaxl = -1e5
  zmmaxl = -1e5
  do i=1,nhalo
     if(rth(i) == RTHmax) then ! get x,y,z coords of most massive halo
        xmmax = posxyz(3,i) 
        ymmax = posxyz(2,i) 
        zmmax = posxyz(1,i) 
        xmmaxl = xmmax
        ymmaxl = ymmax
        zmmaxl = zmmax
        rmmax = sqrt(xmmax**2+ymmax**2+zmmax**2)
        write(*,18) RTHmax,rmmax, xmmax, ymmax, zmmax
18      format(/,'Before rotation largest halo is at',/,&
             'RTH, distance:    ',2(1pe9.3,1x),/,&
             'x, y, z:          ',3(1e10.3,1x),/)

     endif
  enddo

  !rotate so largest halo is in first quadrant to make angles easier
  call mpi_allreduce(xmmaxl,xmmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
  call mpi_allreduce(ymmaxl,ymmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
  call mpi_allreduce(zmmaxl,zmmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
  if(xmmax <0.0)   posxyz(3,:) = -posxyz(3,:)
  if(ymmax <0.0)   posxyz(2,:) = -posxyz(2,:)
  if(zmmax <0.0)   posxyz(1,:) = -posxyz(1,:)
  xmmax = abs(xmmax)
  ymmax = abs(ymmax)
  zmmax = abs(zmmax)

  rmaxtheta = asin(zmmax/sqrt(xmmax**2+zmmax**2)) ! rotate in z-x plane
  xmmaxl      =  xmmax/abs(xmmax)*sqrt(xmmax**2+zmmax**2) 
  rmaxphi   = asin(ymmax/sqrt(ymmax**2+xmmaxl**2))  ! rotate in z-y plane
  
  ! Rotate by angles defined by largest halo
  do i=1,nhalo
     xh =   posxyz(3,i)*cos(rmaxtheta) + posxyz(1,i)*sin(rmaxtheta) 
     zh =  -posxyz(3,i)*sin(rmaxtheta) + posxyz(1,i)*cos(rmaxtheta) 

     yh =  -xh*sin(rmaxphi) + posxyz(2,i)*cos(rmaxphi) 
     xh =   xh*cos(rmaxphi) + posxyz(2,i)*sin(rmaxphi) 
     
     !change (x',0,0) to (0,0,z')
     posxyz(3,i) = zh
     posxyz(2,i) = yh
     posxyz(1,i) = xh 
     
  enddo
  
  ! Double check rotation worked and find distance to largest halo
  ! find x,y,z position and angles of largest halo
  do i=1,nhalo
     if(rth(i) == RTHmax) then
        xmmax = posxyz(3,i) 
        ymmax = posxyz(2,i) 
        zmmax = posxyz(1,i) 
        rmmax = sqrt(xmmax**2+ymmax**2+zmmax**2)
        write(*,19) RTHmax,rmmax, xmmax, ymmax, zmmax
19      format(/,'After rotation largest halo is at',/,&
             'RTH, distance:    ',2(1pe9.3,1x),/,&
             'x, y, z:          ',3(1e10.3,1x),/)
     endif
  enddo


111 continue

  mass = 4. / 3. * mypi * rho_0 * rth ** 3
  
  ! Load table

  call loadmaptable(tablefile)

  ! theta0 and phi0 are the center of fov
  ! here we set the center to lie along z-axis
  theta0 = 0
  phi0 = 0
  theta1 = mypi/2 - theta0
  
  ! Loop over halos

  fsmapl=0
  hpmapl=0
  maxtheta=0
  m=0
  do i=1,nhalo

     xh  = posxyz(3,i)
     yh  = posxyz(2,i) ! Note that x and z are switched
     zh  = posxyz(1,i)           
     vh  = vrad(i)
     
     chih = sqrt(xh**2+yh**2+zh**2)
     redshifth=zofr(chih)

     xh = xh - chihview
     chih = sqrt(xh**2+yh**2+zh**2)

     mh = sqrt(deltacrit(redshifth)/bbps_delta(model)) * mass(i) ! Assuming SIS
     if(mass(i)<mmin) cycle
     m = m+1
     rh  = (3*mh/4/mypi/delta_pksc/rhocrit(redshifth))**(1./3.)
!     call dot(
     call projecthalo_flatsky(xh,yh,zh,rh,chih,mh,redshifth,vh,profile)
     call projecthalo_healpix(xh,yh,zh,rh,chih,mh,redshifth,vh,profile)

  enddo

  nhalol        = m
  maxthetal     = maxtheta
  maxthetamassl = maxthetamass
  maxthetachihl = maxthetachih

  call mpi_reduce(nhalol,nhalo,1,mpi_int,mpi_sum,0,mpi_comm_world,ierr)

  ! Gather maps from processors
  call mpi_allreduce(fsmapl,fsmap,npix**2,mpi_real,mpi_sum,mpi_comm_world,ierr)
  call mpi_allreduce(hpmapl,hpmap,hpnpix, mpi_real,mpi_sum,mpi_comm_world,ierr)

  if(myid==0) then

  write(*,15) nhalo
15 format('Number of halos in map:    ',i0)

  fsmom1 = sum(fsmap**1) /npix / npix
  fsmom2 = sum(fsmap**2) /npix / npix
  fsmom3 = sum(fsmap**3) /npix / npix

  hpmom1 = sum(hpmap**1) / hpnpix
  hpmom2 = sum(hpmap**2) / hpnpix
  hpmom3 = sum(hpmap**3) / hpnpix
  
  write(*,16) fsmom1,fsmom2,fsmom3,hpmom1,hpmom2,hpmom3
16 format('Flatsky moments:    ',3(1pe9.3,1x),/,&
          'Healpix moments:    ',3(1pe9.3,1x))

  write(*,17) sqrt(sum((fsmap-fsmom1)**2)/npix**2),fsmom1,&
              sqrt(sum((hpmap-hpmom1)**2)/hpnpix) ,hpmom1
17 format('Flatsky RMS:                ',1pe9.3,/,&
          'Flatsky Mean:               ',1pe8.2,/,&
          'Healpix RMS:                ',1pe9.3,/,&
          'Healpix Mean:               ',1pe8.2,/)
     

  ! Get dp/dy
  dpdy=0
  do i=0,hpnpix-1
     m = int(hpmap(i,1)/dpdy_dy)+1
     if(m>0.and.m<=ndpdy) dpdy(m)=dpdy(m)+1
  enddo
  dpdy = dpdy / real(hpnpix)

  ! Check mean
  fsmom1=0.0
  do i=1,ndpdy
     fsmom1=fsmom1+(i-0.5)*dpdy_dy*dpdy(i)
  enddo

  do i=2,ndpdy
     dpdy(i)=dpdy(i-1)+dpdy(i)
  enddo  
 
  ilast = indlnb(fileout)

  fileout_fs=fileout(1:ilast)//'_fs.map'
  fileout_hp=fileout(1:ilast)//'_hp.fits'
  fileout_dp=fileout(1:ilast)//'_dp.bin'

  ! P(<y) file
  open(unit=1,file=fileout_dp,form='binary')
  write(1) ndpdy,dpdy_max
  write(1) dpdy
  close(1)

  ! Flat sky binary file
  open(unit=1,file=fileout_fs,form='binary')
  write(1) npix,npix,fov,fov
  write(1) fsmap
  close(1)
     
  ! Healpix FITS file
  hpheader(:)=''
  call add_card(hpheader,'NSIDE',hpnside,'the nside of the map')
  call add_card(hpheader,'ORDERING','RING','the nside of the map')
  call output_map(hpmap,hpheader,fileout_hp)

  endif

  if(myid==0) write(*,*)

  call mpi_finalize(ierr)

  stop

81 format('000',i1)
82 format( '00',i2)
83 format(  '0',i3)
84 format(      i4)

end program pks2map
