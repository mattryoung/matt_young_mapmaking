program hpkvdmain

  use mpivars

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  call hpkvd

  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_finalize(ierr)

contains 

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine hpkvd

  ! ======================================================
  ! DECLARATIONS
  ! ======================================================

  ! INCLUDE MODULES
  use intreal_types
  use mpivars
  use timing_diagnostics
  use input_parameters
  use textlib
  use fftw_interface

  use params
  use ccoordstrat
  use EdeS
  use tabvec
  use seq_name
  use sig
  use evalues
  use etafiles
  use etafiles1
  use etafiles2
  use table
  use particle

  use RandomFieldWrapper
  use SlabToCube

  use memory_management

  use io

  ! PARAMETERS
  real,    parameter :: pi = 4.0 * atan(1.0)
  real,    parameter :: f3p = 4.0/3.0*pi

  ! LATTICE 
  integer  nn(3), nhunt
  real     alattv(3),Sbar(3),Ebar(3,3),Sbar2(3)

  ! FILTERBANK 
  integer  idclus(nclmax)
  real     fcritv(nclmax),Rfclv(nclmax)
  integer  nic, nicz
  ! HALOS
  integer, allocatable :: Npk_eachtask(:), Npk_eachtaskl(:), Npk_begin(:)
  integer, allocatable :: ndump_eachtask(:), ndump_eachtaskl(:), ndump_begin(:)
  
  ! REDSHIFTS
  real, allocatable :: redshifts(:)  

  ! DATE AND TIME 
  character(8)          :: dnt_date
  character(10)         :: dnt_time
  character(5)          :: dnt_zone
  integer, dimension(8) :: dnt_values
  integer initial_time, final_time, count_rate
  integer elapsed_hours, elapsed_minutes, elapsed_seconds
  double precision elapsed_time

  ! MEMORY USAGE
  real grid_mem, fftw_mem, s2c_mem, cat_mem, lat_mem, total_mem
  real local_grid_size, fftw_grid_size

  ! HALO STATISTICS
  real    xpkmin,   xpkmax,   ypkmin,   ypkmax,   zpkmin,   zpkmax
  real    vxpkmin,  vxpkmax,  vypkmin,  vypkmax,  vzpkmin,  vzpkmax
  real    xpkminl,  xpkmaxl,  ypkminl,  ypkmaxl,  zpkminl,  zpkmaxl
  real    vxpkminl, vxpkmaxl, vypkminl, vypkmaxl, vzpkminl, vzpkmaxl
  real    RTHLmin
  real    RTHLminl, RTHLmaxl

  ! RANDOM FIELD PARAMETERS
  character *128 rf_filepktab,rf_outcode
  real           rf_boxsize,rf_ainit
  integer        rf_nmesh,rf_nbuff,rf_ntile,rf_seed

  ! SHADOW BOXES
  logical findpeaks

  double precision D, abar
  

  ! I/O
  character*128 outputfilename
  integer(C_LONG) :: pos_offset
  integer(C_LONG)  :: bytesout
  real(C_FLOAT),  pointer :: outbuffer(:)
  real(C_FLOAT),  target :: RTHLmax
  integer(C_INT), target :: npk_tot
  integer	  :: outnum

  ! ======================================================
  ! EXECUTABLE STATMENTS
  ! ======================================================

  ! SEED IS NOW A COMMANDLINE PARAMETER FOR MC RUNS
  iseedFFT = i4arg(1,13579)

  ! HARD CODES NOT YET INCLUDED IN PARAMETER FILE
  verbose = .false.
  rf_report = 1
  report_memory=.true.

  ! REPORT DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     write(*,11) dnt_values(3),dnt_values(2),&
          dnt_values(1),dnt_values(5),dnt_values(6),&
          dnt_values(7),dnt_zone
     call system_clock(initial_time,count_rate)
  endif

  if(report_memory) &
       call ReportMemory(C_CHAR_"Initial"//C_NULL_CHAR)

  ! READ PARAMETERS AND SET COSMOLOGY
  call read_parameters
  call read_cosmology(Omt,Omnr,Omx,OmB,h,&
       Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)
  call Dlinear_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,&
       fhdmclus,iamcurved,dcurv)

  ! INTIALIZE RANDOM FIELD
  rf_filepktab = pkfile
  rf_outcode   = filein
  rf_boxsize   = dcore_box*(nlx-1) + dL_box 
  rf_nmesh     = n1
  rf_nbuff     = nbuff
  rf_ntile     = nlx
  rf_ainit     = 1.0
  rf_seed      = iseedFFT
  call RandomField_init(rf_filepktab,rf_outcode,rf_boxsize,rf_nmesh,&
       rf_nbuff,rf_ntile,rf_ainit,rf_seed,NonGauss,fNL)

  ! INTIALIZE SLAB TO CUBE
  call SlabToCube_init(nlx,nly,nlz)

  ! MISCELLANEOUS PARAMETERS 
  nn(1)=n1
  nn(2)=n2
  nn(3)=n3
  n1_2=n1/2
  nyq1=n1/2+1
  n1xn2=n1*n2
  nF=nFmax
  alatt   = dL_box /  n1
  boxlen  = alatt * n1
  alatt_2 = (1 / alatt)**2
  ak_min  = 2 * pi / boxlen
  alattv = alatt
  vTHvir0  = 100 * h * sqrt(Omnr)
  vTHvir02 = vTHvir0*vTHvir0
  outnum = 0

  ! READ FILTERBANK AND MAKE BOXES
  call read_filterbank(nic,idclus,fcritv,Rfclv,Rfclmax)

  ! INITIALIZE LATTICE HUNT ARRAY
  ! nhunt = Rfclmax / alatt 
  !GS 23/02/2016 - see ir2min comment                                           
  nhunt = 48.0/alatt
  call icloud(npartmax,nhunt)

  ! REPORT MEMORY BUDGET
  fftw_grid_size  = total_local_sizes
  local_grid_size = (n1+2)*n1**2

  fftw_mem = 16.*fftw_grid_size/1024**3     ! slab(4)
  s2c_mem  = 4. *local_nz*n1**2/1024**3     ! sendbuff(4)
  grid_mem = 21.*local_grid_size/1024**3    ! rho(4),eta(12),mask(1),recvbuff(4)
  cat_mem  = 92.*(2*npkmax+npkmaxl)/1024**3 ! catalogs(184)
  lat_mem  = 74.*npartmax/1024**3           ! nshell(4),Sshell(12),SRshell(36),
                                            ! Fbar(4),rad(4),indx(4),ixsvec(6)
                                            ! irs2(4)
                                            ! --------------------------------
                                            ! 
  total_mem   = fftw_mem + grid_mem + &     ! FFTW + GRID + CATALOG + LATTICE +
                cat_mem + lat_mem + &       ! S2C BUFFER
                s2c_mem

  if(myid==0) then
        write(*,*)
        write(*,*) '  FFTW memory:        ',fftw_mem,   ' GB'
        write(*,*) '  Grid memory:        ',grid_mem,   ' GB'
        write(*,*) '  S2C memory:         ',s2c_mem,    ' GB'
        write(*,*) '  Catalog memory:     ',cat_mem,    ' GB'
        write(*,*) '  Lattice hunt memory:',lat_mem,    ' GB'
        write(*,*) '  Total memory:       ',total_mem  ,' GB'
	write(*,*) 
  endif

  ! ALLOCATE MASK AND FFT ARRAYS
  allocate(mask(n1,n2,n3))

  ! LOCAL
  call fftw_allocate_serial( delta,    deltac,   deltap)
  call fftw_allocate_serial(delta_u, deltac_u, delta_up) ! to save unsmoothed
  call fftw_allocate_serial(   etax,    etaxc,    etaxp)
  call fftw_allocate_serial(   etay,    etayc,    etayp)
  call fftw_allocate_serial(   etaz,    etazc,    etazp)

  ! GLOBAL
  call fftw_allocate_parallel(deltag, deltagc, deltagp)
  call fftw_allocate_parallel( etaxg,  etaxgc,  etaxgp)
  call fftw_allocate_parallel( etayg,  etaygc,  etaygp)
  call fftw_allocate_parallel( etazg,  etazgc,  etazgp)

  if(report_memory) &
       call ReportMemory(C_CHAR_"After grid allocations"//C_NULL_CHAR)

  ! RANDOM FIELD REALIZATION
  timing_diagnostics_code='field realization'
  call timer_begin

  call RandomField_make(-1, deltag, deltagc, deltag)
  if(ioutfield>0) call RandomField_Output
  call RandomField_make(-2,  etaxg,  etaxgc, deltag) !third field is the 
  if(ioutfield==2) call RandomField_Output
  call RandomField_make(-3,  etayg,  etaygc, deltag) !one to convolve with
  if(ioutfield==2) call RandomField_Output  
  call RandomField_make(-4,  etazg,  etazgc, deltag) 
  if(ioutfield==2) call RandomField_Output
  call timer_end  

  ! ALLOCATE REDSHIFTS
  allocate(redshifts(num_redshifts))  

  ! MAKE HOMOGENEOUS ELLIPSOID TABLES
  timing_diagnostics_code='Collapse tables'
  call timer_begin
  call make_F_zvir_table
  call inithom_ellipse_params
  call makehom_ellipse_tab
  if(myid==0) call timer_end

  ! MAKE BOXES AND BOX TO PROCESS MAPPING
  m = 0
  call allocate_boxes
  call make_boxes
  do i=1,nboxes

     ! HERE THEY ARE ABSOLUTE VALUES OF BOX NEAR CORNERS
     xbx = abs(xbox(i)) - dcore_box / 2
     ybx = abs(ybox(i)) - dcore_box / 2
     zbx = abs(zbox(i)) - dcore_box / 2

     rbx = sqrt(xbx**2+ybx**2+zbx**2)
     abx = afn(rbx)
     zbx = 1 / abx - 1

     if(zbx > maximum_redshift.and.ievol==1) cycle

     m = m + 1
     boxlist(m) = i
     
  enddo
 
  nboxes = m
  neach  = int(nboxes/ntasks)+1

  ! GET GRID DENSITY ASSIGNMENT KERNEL TABLE
  call atab4

  ! INITIALIZE REDSHIFTS

  dz = (maximum_redshift-global_redshift)/num_redshifts
  do i=1,num_redshifts
    redshifts(i) = global_redshift + (i-1)*dz
  enddo
  if(num_redshifts>1) ievol=0

  ! INITIALIZE NUMBER OF PEAKS 
  Npk  = 0
  ndump=0

  ! ALLOCATE HALOS
  call allocate_halos(num_redshifts)  

  ! ======================================================
  ! BEGIN MAIN PEAK PATCH CALCULATION
  ! ======================================================
  ! LOOP OVER BOXES
  do ibox=myid+1,ntasks*neach,ntasks

     findpeaks = .true.
     if(ibox>nboxes) then
       jbox = boxlist(ibox-nboxes)
       findpeaks = .false.
     else
       jbox = boxlist(ibox)
     endif

     xbx = xbox(jbox)
     ybx = ybox(jbox)
     zbx = zbox(jbox)

     rbx = sqrt(xbx**2+ybx**2+zbx**2)
     abx = afn(rbx)
     rdbx = 1/abx - 1

     delta = 0
     etax  = 0
     etay  = 0 
     etaz  = 0

     ! SLAB TO CUBE EXCHANGE

     call timer_begin
     timing_diagnostics_code='Slab to cube exchange'

     call get_cube_indices(xbx,ybx,zbx,nx,ny,nz)
     call SlabToCube_map(nx,ny,nz,nlx,nly,nlz)

     call SlabToCube_exchange(delta, deltag)
     call SlabToCube_exchange( etax,  etaxg)
     call SlabToCube_exchange( etay,  etayg)
     call SlabToCube_exchange( etaz,  etazg)

     call timer_end  

     if(findpeaks) then

     ! LOOP OVER REDSHIFTS
     do ired=1,num_redshifts

     mask = 0

     ! MULTIPLY BY LINEAR GROWTH FACTOR

     if(ievol==1) then
      cen  = 0.5*(n1+1)
      D = 0
      abar = 0
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(alatt,cen,xbx,ybx,zbx,delta,etax,etay,etaz) &
!$OMP REDUCTION(+:abar,D)
        do k=1,n3
           zc = zbx + alatt * (k - cen)
           do j=1,n2
              yc = ybx + alatt * (j - cen)
              do i=1,n1
                 xc = xbx + alatt * (i - cen)
                 rc = sqrt(xc**2+yc**2+zc**2)

                 a  = afn(rc)
                 growth_factor = Dlinear(a,rc,HD_a,D_a)

		 delta(i,j,k) = delta(i,j,k) * growth_factor
		 etax(i,j,k) = etax(i,j,k) * growth_factor
		 etay(i,j,k) = etay(i,j,k) * growth_factor
		 etaz(i,j,k) = etaz(i,j,k) * growth_factor

                 abar = abar + a
	         D    = D    + growth_factor

              enddo
           enddo
        enddo

        abar = abar / n1 / n2 / n3
	D    = D    / n1 / n2 / n3

        z    = 1 / abar - 1

        !GS 23/02/2016 add empirically derived redshift evolving filterbank
        nicz = min( nic, int(nic / (1+z)**(1./2) + 3.5) )
     else
        z   = redshifts(ired)
        a   = 1 / (1+z)
        rc  = 0.0
        D   = Dlinear(a,rc,HD_a,D_a)

        delta = delta * D
        etax = etax * D
        etay = etay * D
        etaz = etaz * D

        nicz = nic
     endif

     ! SET TABLES FOR THIS REDSHIFT
     call inithom_ellipse_redshift(1.0)

     if(myid==0) write(*,101) ibox,nboxes,D,z,Npk

     timing_diagnostics_code='Smoothing'
     call timer_begin

     !ADDED BY G.S 16/07/2015
     !SMOOTH WITH "MINIMAL" FILTER R_G = 0.5 * alatt as in BM1 Appendix E
     !TO DAMP HIGH-FREQUENCY NOISE 
     !MULT R_G BY 2 TO ACCOUNT FOR THE DIFFERENCE IN TOP-HAT TO GAUSSIAN
     !FILTERING ADDED IN SMOOTH_FIELD

     call fftwf_execute_dft_r2c(plan_s, delta, deltac_u)
     !call smooth_field(nyq1,ak_min,n1_2,0.5*alatt*2,fknyq,speq)
     !deltac_u = deltac
     ! SAVE "UNSMOOTHED" (MINIMALLY SMOOTHED) DENSITY IN RECIPROCAL SPACE

     ! LOOP OVER SMOOTHING SCALES
     if(myid==0) write(*,102) 
     Nprev=Npk
     nbad=0

     do ic=nicz,1,-1

        !SMOOTHING IS NOW IN FUNCTION
        call smooth_field(nyq1,ak_min,n1_2,Rfclv(ic),fknyq,speq)

        ! FFT TO REAL SPACE
        call fftwf_execute_dft_c2r(iplan_s, deltac, delta)
        delta = delta / real(n1)**3
        sigrho=sqrt(sum(delta(1:n1,:,:)**2)/n1/n2/n3)
        if(myid==0) write(*,103) Rfclv(ic),sigrho
!        write(*,105) maxval(delta), Rfclv(ic)
!        write(*,*) maxval(delta), Rfclv(ic)

        f_v = fsc_tabfn_of_ZZ(1+z)
	! FIND PEAKS
        Npk_in=Npk
        call get_pks(f_v,Npk,xbx,ybx,zbx,alattv,ired, Rfclv(ic))
        if (Npk.ge.Npkmaxl) then
           if(myid==0)  write(*,*) 'Reached cluster limit ',Npkmaxl
           call mpi_finalize(ierr)
           stop
        endif
        Npk_new=Npk-Npk_in
        if(myid==0) write(*,104) Npk_new
!        write(*,104) Npk_new
        if(Npk_new.lt.1) cycle

        do jpp=Npk_in+1,Npk
           iwas(jpp,ired)=ic
        enddo
          
     enddo
     if(myid==0) call timer_end
     
     ! END LOOP OVER SMOOTHING SCALES
     if(myid==0) write(*,107) Npk-Nprev

     ! SET DELTA BACK TO UNSMOOTHED VALUE
     call fftwf_execute_dft_c2r(iplan_s, deltac_u, delta)
     delta = delta / real(n1)**3

     ! LOOP OVER PEAKS
     timing_diagnostics_code='Ellipsoidal collapse'
     call timer_begin
     do jpp=Nprev+1,Npk
        if(mod((jpp-Nprev+1),5000)==0.and.myid==0) then
           write(*,108) jpp-Nprev+1
        endif
        ic=iwas(jpp,ired)
        jp=lagrange(jpp,ired)

        if (ievol.eq.1) then
           radpk=sqrt(xpk(jpp,ired)**2+ypk(jpp,ired)**2+zpk(jpp,ired)**2)
           chipk=radpk
           apk=afn(chipk)
           redshiftpk=1/apk - 1
        else
           redshiftpk=z
        endif
        redshiftpk=0

        Fnupk=fsc_tabfn_of_ZZ(1+redshiftpk)
        Fevpk=0.
        Fpvpk=0.
        ir2min=min( int((1.75*Rfclv(ic)/alatt)**2), int((48.0/alatt-1)**2)) !GS hardcode 45Mpc instead of 
        if(rmax2rs .gt. 0.0) then                                           !Rfclv(nic) so don't 
           npart = int(f3p*(rmax2rs*Rfclv(ic)/alatt)**3)                    !waste smoothing
           npart = min(npart,npartmax)
        else
           npart = npartmax
        endif

        ! HOMOGENEOUS ELLIPSOID CALCULATION
        call get_homel(npart,jp,alatt,ir2min,&
             Sbar,RTHL,Srb,1+redshiftpk,Fnupk,Fevpk,&
             Fpvpk,Ebar,Sbar2,Rfclv(ic))
        report_memory=.false.
        RTHLv(jpp,ired)=RTHL*alatt
        if(RTHLv(jpp,ired)<=0) then
           RTHLv(jpp,ired)=0.0
           vTHvir(jpp,ired)=0.0
           iwas(jpp,ired)=-10000
           ndump=ndump+1
        else
           vE2 = vTHvir02 * Fnupk * Srb * RTHLv(jpp,ired)**2
           vTHvir(jpp,ired)=sqrt(abs(vE2))
           if(vE2.lt.0.0) vTHvir(jpp,ired)=-vTHvir(jpp,ired)              

           ! PUT SHEAR EIGENVALUES F_nu,F_ev,F_pv WHERE
           !
           !     lam1 = F_nu/3( 1 + 3*F_ev + F_pv )
           !     lam2 = F_nu/3( 1 - 2*F_pv )
           !     lam3 = F_nu/3( 1 - 3*F_ev + F_pv )
           !

           Fcollv(jpp,ired)=Fnupk
           F_ev(jpp,ired)=Fevpk
           F_pv(jpp,ired)=Fpvpk

           if(ioutshear.eq.1) then
	      strain_bar(1,jpp,ired)=Sbar2(1)
	      strain_bar(2,jpp,ired)=Sbar2(2)
	      strain_bar(3,jpp,ired)=Sbar2(3)

	      !c PB 2015.06.18 - Switched to local shear measurement
              !c strain_bar(1,jpp,ired)=Ebar(1,1)
              !c strain_bar(2,jpp,ired)=Ebar(2,2)
              !c strain_bar(3,jpp,ired)=Ebar(3,3)
              !c strain_bar(4,jpp,ired)=Ebar(2,3)
              !c strain_bar(5,jpp,ired)=Ebar(1,3)
              !c strain_bar(6,jpp,ired)=Ebar(1,2)
           endif

           vxpk(jpp,ired)=Sbar(1)
           vypk(jpp,ired)=Sbar(2)
           vzpk(jpp,ired)=Sbar(3)

        endif
        
        if (iwas(jpp,ired).gt.-10000) iwas(jpp,ired)=idclus(ic)
        if (idoboxf.eq.1) lagrange(jpp,ired)=jbox

     enddo
     if(myid==0) call timer_end

     ! END LOOP OVER PEAKS
  
     if(ievol==0) then 
        delta = delta / D
        etax  = etax  / D
        etay  = etay  / D
        etaz  = etaz  / D
     endif

  enddo
  ! END LOOP OVER REDSHIFTS

  endif

  enddo
  ! END LOOP OVER BOXES
  ! ======================================================
  ! END MAIN PEAK PATCH CALCULATION
  ! ======================================================

  call mpi_barrier(mpi_comm_world,ierr)

  report_memory=.true.

  ! FREE ALL GRID MEMORY

  if(report_memory) &
       call ReportMemory(C_CHAR_"Before grid deallocations"//C_NULL_CHAR)

  ! LOCAL 
  call fftwf_free(deltap)
  call fftwf_free(etaxp)
  call fftwf_free(etayp)
  call fftwf_free(etazp)
  
  ! GLOBAL
  call fftwf_free(deltagp)
  call fftwf_free(etaxgp)
  call fftwf_free(etaygp)
  call fftwf_free(etazgp)

  if(report_memory) &
       call ReportMemory(C_CHAR_"After grid deallocations"//C_NULL_CHAR)

  ! FIND OFFSET AND TOTAL NUMBER OF PEAKS
  allocate(Npk_eachtaskl(0:ntasks-1))
  allocate(Npk_eachtask (0:ntasks-1))
  allocate(Npk_begin    (0:ntasks-1))
  allocate(ndump_eachtaskl(0:ntasks-1))
  allocate(ndump_eachtask (0:ntasks-1))
  allocate(ndump_begin    (0:ntasks-1))

  Npk_begin           = 0
  Npk_eachtaskl       = 0
  Npk_eachtaskl(myid) = Npk-ndump

  npk_begin             = 0
  ndump_eachtaskl       = 0
  ndump_eachtaskl(myid) = ndump

  call mpi_allreduce(Npk_eachtaskl,Npk_eachtask,ntasks,mpi_integer,mpi_sum,&
                     mpi_comm_world,ierr)
  call mpi_allreduce(ndump_eachtaskl,ndump_eachtask,ntasks,mpi_integer,mpi_sum,&
                     mpi_comm_world,ierr)

  Npk_begin(0) = 0
  ndump_begin(0) = 0
  do i=1,ntasks-1
    Npk_begin(i)=Npk_eachtask(i-1)+Npk_begin(i-1)
    ndump_begin(i)=ndump_eachtask(i-1)+ndump_begin(i-1)
  enddo

  !PB 2015.06.18
  if(ioutshear.eq.0) then
      outnum = 7
  endif
  if(ioutshear.eq.1) then
      outnum = 10
  endif 
  

  Npk_tot    = Npk_begin(ntasks-1) + Npk_eachtask(ntasks-1)
  ndump_tot    = ndump_begin(ntasks-1) + ndump_eachtask(ntasks-1)

  ! LOAD OUTPUT BUFFERS
  fileout=trim(fileout)//'.'//str(iseedFFT)

  if(report_memory) &
       call ReportMemory(C_CHAR_"Before output buffer allocation"//C_NULL_CHAR)

  outbufferp = fftwf_alloc_real(int(outnum*(Npk-ndump), C_SIZE_T))
  call c_f_pointer(outbufferp,outbuffer,[outnum*(Npk-ndump)])

  if(report_memory) &
       call ReportMemory(C_CHAR_"After output buffer allocation"//C_NULL_CHAR)

  m=0
  do ired=1,num_redshifts
     do i=1,Npk
        !G.S 24/09/2015 only write out good peaks (Npk-ndump)
        if (RTHLv(i,ired) > 0.0) then  
           outbuffer(m+1) =   xpk(i,ired)
           outbuffer(m+2) =   ypk(i,ired)
           outbuffer(m+3) =   zpk(i,ired)
           outbuffer(m+4) =  vxpk(i,ired)
           outbuffer(m+5) =  vypk(i,ired)
           outbuffer(m+6) =  vzpk(i,ired)
           outbuffer(m+7) = RTHLv(i,ired)
        
           if(ioutshear.eq.1) then
              outbuffer(m+8) = strain_bar(1,i,ired)
              outbuffer(m+9) = strain_bar(2,i,ired)        
              outbuffer(m+10) = strain_bar(3,i,ired)        
           endif
           
           m = m + outnum
        endif
     enddo
  enddo

  ! GATHER PEAK INFORMATION
  xpkminl   = minval(xpk)
  xpkmaxl   = maxval(xpk)
  ypkminl   = minval(ypk)
  ypkmaxl   = maxval(ypk)
  zpkminl   = minval(zpk)
  zpkmaxl   = maxval(zpk)

  vxpkminl  = minval(vxpk)
  vxpkmaxl  = maxval(vxpk)
  vypkminl  = minval(vypk)
  vypkmaxl  = maxval(vypk)
  vzpkminl  = minval(vzpk)
  vzpkmaxl  = maxval(vzpk)

  RTHLminl = minval(RTHLv)
  RTHLmaxl = maxval(RTHLv)

  call mpi_reduce(xpkminl,xpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(xpkmaxl,xpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(ypkminl,ypkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(ypkmaxl,ypkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(zpkminl,zpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(zpkmaxl,zpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  call mpi_reduce(vxpkminl,vxpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vxpkmaxl,vxpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(vypkminl,vypkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vypkmaxl,vypkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(vzpkminl,vzpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vzpkmaxl,vzpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  call mpi_reduce(RTHLminl,RTHLmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(RTHLmaxl,RTHLmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  ! WRITE HEADER IF myid == 0
  bytesout   = 4
  if(myid==0) then

    ! WRITE npk_tot
    pos_offset = 0
    call myfwrite(trim(fileout)//c_null_char,&
                  c_loc(npk_tot), bytesout, pos_offset)  

    ! WRITE RTHLmax
    pos_offset = pos_offset + bytesout
    call myfwrite(trim(fileout)//c_null_char,&
                  c_loc(RTHLmax), bytesout, pos_offset)

  endif 

  ! WRITE PEAKS
  pos_offset = int(4,8)*(2+outnum*(Npk_begin(myid))) 
  bytesout   = 4*outnum*(Npk-ndump)
  if(Npk  >0) call myfwrite(trim(fileout)//c_null_char,&
                            outbufferp, bytesout, pos_offset)

  ! REPORT PEAK INFORMATION
  if(myid==0) then
     write(*,301) Npk_tot
     write(*,302) ndump_tot
     write(*,303) RTHLmin,RTHLmax
     write(*,304) xpkmin,xpkmax,ypkmin,ypkmax,zpkmin,zpkmax
     write(*,305) vxpkmin,vxpkmax,vypkmin,vypkmax,vzpkmin,vzpkmax
  endif
  
  ! DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     call system_clock(final_time,count_rate)

     elapsed_time    = (final_time - initial_time) / count_rate
     elapsed_hours   = int(elapsed_time / 3600)
     elapsed_minutes = int((elapsed_time - elapsed_hours * 3600)/60)
     elapsed_seconds = int(elapsed_time - elapsed_hours * 3600 &
                                        - elapsed_minutes * 60)
     
     write(*,401) dnt_values(3),dnt_values(2),dnt_values(1),&
       dnt_values(5),dnt_values(6),dnt_values(7),&
       dnt_zone,elapsed_hours,elapsed_minutes,elapsed_seconds

  endif

  return

 11 format(/,3x,61('-'),/,3x,'Peak-patch HPKVD running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)
101 format(/,3x,'Box ',i0,' of ',i0,/,3x,&
         'D = ',f5.3,/,3x,'z = ',f5.2,/,3x,'Npks = ',i0)
102 format(/,5x,'Smoothing:')
103 format(  5x,'R = ',f7.3,' Mpc: sigma = ',f5.3)            
104 format(  10x,i0, ' peaks found')            
105 format(  10x, ' max delta(R), R = ',f7.3,f7.3)
107 format(/,5x,'Homogeneous Ellipsoid running on ',i0,&
         ' peaks:')
108 format(  5x,'Peak',1x,i0,1x,'done')
201 format(3x,/, 'Gathering peak information on process 0')
301 format(/,3x, 'Total number of peaks kept:   ',11x,i0)
302 format(3x, 'Total number of peaks dumped: ',11x,i0)
303 format(/,3x,   'Minimum and maximum of peak radii:',&
           /,5x,1('(',1pe13.6',',1pe13.6,')',2x,/))
304 format(3x,   'Minimum and maximum of coordinates:',&
           /,5x,3('(',1pe13.6',',1pe13.6,')',/,5x))
305 format(3x,   'Minimum and maximum of displacements:',&
           /,5x,3('(',1pe13.6',',1pe13.6,')',/,5x))

401 format(/,3x,61('-'),/,3x,&
         'Peak-patch HPKVD finished running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         40x,'Total run time: ',i2.2,':',i2.2,':',i2.2,/&
         3x,61('-'),/)
                                    
end subroutine hpkvd


! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!SMOOTHES COMPLEX FIELD ON SCALE Rfclv AND RETURNS        

subroutine smooth_field(nyq1,ak_min,n1_2,Rfclv,fknyq,speq)

  use input_parameters
  use mpivars
  use arrays
  use fftw_interface

  sigrho=0.0
  fknyq=(nyq1-1)*ak_min
!$OMP PARALLEL DO &                                                  
!$OMP DEFAULT(FIRSTPRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(nyq1,ak_min,n1_2,Rfclv,delta,deltac,fknyq,speq)
  do k=1,n3
     kk=k
     if(k.ge.nyq1) kk=n3+2-k
     fkz=(kk-1)*ak_min

     do j=1,n2
        jj=j
        if(j.ge.nyq1) jj=n2+2-j
        fky=(jj-1)*ak_min

        do i=1,n1_2+1
           ii=i
           if(i>n1_2) ii=n1+2-i
           fkx=(ii-1)*ak_min

           ! CONVOLVE WITH GAUSSIAN OR TOP-HAT FILTER                                                
           fkR=sqrt(fkx**2+fky**2+fkz**2)*Rfclv
           if(wsmooth==0) then
              fkR = fkR/2
              w = wgauss(fkR)
           elseif(wsmooth==1) then
              w = wth(fkR)
           endif
           if(fkR==0) w=1
           deltac(i,j,k) = deltac_u(i,j,k) * w

        enddo

     enddo
  enddo

  return
end subroutine smooth_field


! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine get_cube_indices(xbx,ybx,zbx,i,j,k)

  USE intreal_types
  use input_parameters
  use arrays
  use mpivars
  
  real xmin,ymin,zmin

  xmin = cenx - nlx/2.*dcore_box
  ymin = ceny - nly/2.*dcore_box
  zmin = cenz - nlz/2.*dcore_box

  i = int((xbx - xmin)/dcore_box)+1
  j = int((ybx - ymin)/dcore_box)+1
  k = int((zbx - zmin)/dcore_box)+1

  return

end subroutine get_cube_indices

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

real function wgauss(x)

  implicit none
  real x
    
  wgauss = exp(-(x**2)/2.)

  return 

end function wgauss

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

real function wth(x)

  implicit none
  real x

  wth = 3.*(sin(x)-x*cos(x))/x**3

  return 

end function wth

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine read_external_field(xbx,ybx,zbx,&
     ibox,speq,speqx,speqy,speqz)

  USE intreal_types
  use input_parameters
  use arrays
  use mpivars
  
  character*128 filename,filebox
  
  integer ibox
  real    xbx,ybx,zbx
  
  complex speq(  :,:)
  
  complex speqx(  :,:),speqy(  :,:),speqz(  :,:)

  integer nxc,nyc,nzc,ntot
  real    boxsize, boxext
  
  ! FIND WHERE THIS BOX IS IN INPUT FIELD
  ntot    = nlx * ( n1 - 2 * nbuff )
  boxsize = dcore_box * nlx 
  boxext  = boxsize * next / ntot
  dx      = boxext / next 
  deltax  = xbx + boxext / 2
  deltay  = ybx + boxext / 2
  deltaz  = zbx + boxext / 2
  nxc     = int( deltax / dx + 1.1)
  nyc     = int( deltay / dx + 1.1)
  nzc     = int( deltaz / dx + 1.1)
  
  ! READ EXTERNAL FIELD 
  
  filename='Fvec_'//filein
  F = delta
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  filename='etax_'//filein
  F = etax
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  filename='etay_'//filein
  F = etay
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  filename='etaz_'//filein
  F = etaz
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  ! eta = -displacement, so need to multiply by -1
  etax = -etax
  etay = -etay
  etaz = -etaz
  
  return
  
end subroutine read_external_field

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine readsubbox(filename,nxc,nyc,nzc,ibox)
  
  use mpivars
  use omp_lib

  use input_parameters
  use arrays
  use endian_utility

  character*128 filename

  integer nxc,nyc,nzc
  integer nx1,nx2,ny1,ny2,nz1,nz2
  integer ii,jj,kk
  real dum
  integer*8 m,ioffset,joffset,koffset,ne
  integer, dimension(3) :: sizes,sizese,start
  integer parallel_read

  ne=next
  ns=n1

  nx1 = nxc - ns/2 
  nx2 = nxc + ns/2 - 1
  ny1 = nyc - ns/2 
  ny2 = nyc + ns/2 - 1 
  nz1 = nzc - ns/2 
  nz2 = nzc + ns/2 - 1 

  ! READ INPUT DATA
  parallel_read = 0
  if(parallel_read == 0) then
     
     offset = int(4,8)*(ibox-1)*n1**3+1
     open(unit=30,file=filename,access='stream')
     read(30,pos=int(offset)) F
     close(30)
     
  elseif(parallel_read == 1) then

     sizese=ne
     sizes=n1
     start(1)=nx1-1
     start(2)=ny1-1
     start(3)=nz1-1     

     call mpi_type_create_subarray(3,sizese,sizes,start,&
          mpi_order_fortran,mpi_real,newtype,ierr)
     call mpi_type_commit(newtype,ierr)
     call mpi_file_open(mpi_comm_world,filename,mpi_mode_rdonly,&
          mpi_info_null,ifile,ierr)  
     call mpi_file_set_view(ifile,0,mpi_real,newtype,"native",&
          mpi_info_null,ierr)
     call mpi_file_read_all(ifile,F,n1*n2*n3,mpi_real,&
          mpi_status_ignore,ierr)
     call mpi_file_close(ifile,ierr)

  else

     offset = int(4,8)*(ibox-1)*n1**3
     call mpi_file_open(mpi_comm_world,filename,mpi_mode_rdonly,&
          mpi_info_null,ifile,ierr)  
     call mpi_file_set_view(ifile,offset,mpi_real,newtype,"native",&
          mpi_info_null,ierr)
     call mpi_file_read(ifile,F,n1*n2*n3,mpi_real,&
          mpi_status_ignore,ierr)
     call mpi_file_close(ifile,ierr)

  endif

  if(big_endian()) then ! ASSUMES INPUT DATA IS LITTLE ENDIAN
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(F) 
    do k=1,n3
      do j=1,n2
        do i=1,n1
	  F(i,j,k)=swap_endian(F(i,j,k))
	enddo
      enddo
    enddo
  endif

  return

end subroutine readsubbox
 
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine make_boxes

 use mpivars
 use arrays
 USE intreal_types
 use input_parameters
 use mpivars
 use seq_name

 nboxes=nlx*nly*nlz

 i=0
 cen1=0.5*float(nlx+1)
 cen2=0.5*float(nly+1)
 cen3=0.5*float(nlz+1)
 do k3=1,nlz
    zBh=cenz+(k3-cen3)*dcore_box
    do k2=1,nly
       yBh=ceny+(k2-cen2)*dcore_box
       do k1=1,nlx
          xBh=cenx+(k1-cen1)*dcore_box

          i=i+1
          xbox(i)=xBh
          ybox(i)=yBh
          zbox(i)=zBh
          
       enddo
    enddo
 enddo

 return

end subroutine make_boxes

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine read_filterbank(nic,idclus,fcritv,Rfclv,Rfclmax)
  USE intreal_types
  use input_parameters

  dimension fcritv(*),Rfclv(*),idclus(*)

  open(1,file=filterfile,status='old',form='formatted')

  read(1,*) nic

  Rfclmax=0.0
  do ic=1,nic
     read(1,*) idclus(ic),fcritv(ic),Rfclv(ic)
     Rfclmax=max(Rfclmax,Rfclv(ic))
  enddo
  close(1)

  return

end subroutine read_filterbank

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine read_parameters

  use input_parameters
  use arrays

  ! READ INPUT PARAMETERS

  open(unit=1,file="hpkvd_params.txt")
  read(1,*) inputdens,idovoid,ioutshear,&
       global_redshift,maximum_redshift,num_redshifts,&
       Omx,OmB,Omhdm,fhdmclus,Omvac,Omcurv,h,gmfmin,gmfmax,&
       ifil,&
       rfmin,rfmax,&
       lkgaus,SFTon,netaf,FFTon,&
       gmLBFFT,gmUBFFT,&
       (nkradv(i),nkangv(i),lsamlogv(i),ensampv(i),alphv(i),&
       gmLBv(i),gmUBv(i),iseedfv(i),i=1,2),&
       iwantdkFFT,iseed_dkfft,iwrap,&
       filter_RF,&
       igetETA,&
       bZfin,bZinit,&
       ivol,nlx,nly,nlz,&
       dcore_box,dL_box,cenx,ceny,cenz,&
       nbuff,next,&
       ievol,igetD2F,&
       f_camp,&
       maskprev,igethom_strat,iwant_evpvtabhpkvd,Nfsc,ifclinv,&
       Fsmin,Fsmax,zinit_fac_ell,devtab,dpv_evtab,evtabmax,&
       pv_evtabmax,&
       ichange_pspec,wsmooth,rmax2rs,ioutfield, NonGauss, fNL
  read(1,'(a)') filein
  read(1,'(a)') pkfile
  read(1,'(a)') filterfile
  read(1,'(a)') fileout
  read(1,'(a)') etabfile
  close(1)

  ! WRITE INPUT PARAMETERS IF VERBOSE=TRUE

  if(verbose.and.myid==0) &
  write(*,*) inputdens,idovoid,ioutshear,&
       global_redshift,maximum_redshift,&
       Omx,OmB,Omhdm,fhdmclus,Omvac,Omcurv,h,gmfmin,gmfmax,&
       ifil,&
       rfmin,rfmax,&
       lkgaus,SFTon,netaf,FFTon,&
       gmLBFFT,gmUBFFT,&
       iseedFFT,&
       (nkradv(i),nkangv(i),lsamlogv(i),ensampv(i),alphv(i),&
       gmLBv(i),gmUBv(i),iseedfv(i),i=1,2),&
       iwantdkFFT,iseed_dkfft,iwrap,&
       filter_RF,&
       igetETA,&
       bZfin,bZinit,&
       ivol,nlx,nly,nlz,&
       dcore_box,dL_box,cenx,ceny,cenz,&
       nbuff,next,&
       ievol,igetD2F,&
       f_camp,&
       maskprev,igethom_strat,iwant_evpvtabhpkvd,Nfsc,ifclinv,&
       Fsmin,Fsmax,zinit_fac_ell,devtab,dpv_evtab,evtabmax,&
       pv_evtabmax,ichange_pspec,wsmooth,rmax2rs,ioutfield
  if(verbose.and.myid==0) then
     write(*,'(a)') trim(filein)
     write(*,'(a)') trim(pkfile)
     write(*,'(a)') trim(filterfile)
     write(*,'(a)') trim(fileout)
     write(*,'(a)') trim(etabfile)
  endif

end subroutine read_parameters

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

end program hpkvdmain
