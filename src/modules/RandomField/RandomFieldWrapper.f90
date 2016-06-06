module RandomFieldWrapper

  !=======================================================================
  !
  ! GENERATES UNIGRID INITIAL CONDITIONS FOR COSMOLOGICAL SIMULATIONS
  !
  !                                                AUTHOR: MARCELO ALVAREZ
  !                                                LAST EDIT:     02.10.15
  !
  !=======================================================================

contains

  subroutine RandomField_init(filepktab_in,outcode_in,boxsize_in,nmesh_in,&
       nbuff_in,ntile_in,ainit_in,seed_in,NonGauss_in,fNL_in)

    !=====================================================================
    ! DECLARATIONS BEGIN  
    
    ! INCLUDE NECESSARY MODULES
    use grid
    use tiles
    use mpivars
    use gaussian_field
    use pktable
    use globalvars
    use fftw_interface

    ! IMPLICIT STATEMENT    
    implicit none

    double precision, parameter :: pi=3.1415926535897932384626433832795_dp

    ! INPUT PARAMETERS
    character *128 filepktab_in, outcode_in
    real           boxsize_in, ainit_in, fNL_in
    integer        nmesh_in, nbuff_in, ntile_in, seed_in, NonGauss_in

    ! OUTPUT PARAMETERS
    integer total_local_sizes_out
    integer fieldid
    
    ! DECLARATIONS END
    !=====================================================================

    !=====================================================================
    ! EXCUTABLE BEGIN

    filepktab = filepktab_in
    outcode   = outcode_in
    boxsize   = boxsize_in
    nmesh     = nmesh_in
    nbuff     = nbuff_in
    ntile     = ntile_in
    ainit     = ainit_in
    seed      = seed_in
    NonGauss  = NonGauss_in
    fNL       = fNL_in
    nsub = nmesh - 2 * nbuff
    n    = nsub * ntile + 2 * nbuff    

    if(mod(n,ntasks).ne.0) then
       if(myid==0) then
         write(*,*) '  ERROR: Number n not evenly divisible by ntasks in FFTW'
         write(*,*) '         exiting ...'
       endif
       call mpi_finalize(ierr)
       stop
    endif

    z = 1./ainit - 1.

    n12  = n / 2
    n21  = n12 + 1
    n2p1 = 2 * (n12 + 1)
    
    dk   = 2 * pi / boxsize
    kmax = n * dk
    d3k  = dk**3

    ! INITIALIZE FFTW  
    nfft   = n 
    nfft_s = nmesh
    call fftw_initialize()

    ! READ POWER SPECTRUM FROM TABLE
    call read_pktable
    pksav    = sqrt(pksav_in*d3k*n**3)
    tsav     = sqrt(tsav_in*d3k*n**3)
    pkchisav = sqrt(pkchisav_in*d3k*n**3)
    return

  end subroutine RandomField_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine RandomField_make(fieldid,fieldp,fieldpc,fieldpg)

    !=====================================================================
    ! DECLARATIONS BEGIN  
    
    !---------------------------------------------------------------------
    ! INCLUDE NECESSARY MODULES
    !---------------------------------------------------------------------
  
    use grid
    use mpivars
    use gaussian_field
    use pktable
    use globalvars
    !for spherical profile
    use fftw_interface
    use mpivars
    use openmpvars
    use tiles
    use chi2zeta

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------
    
    implicit none
    
    integer fieldid
    !for single spherical overdensity
    real alpha
    real ro, rmax, ar, dr, ff, rx,ry,rz, fcrit
    real*8 lavg,avg,lsig,sig

    real(C_FLOAT),            pointer :: fieldp(:,:,:)
    real(C_FLOAT),            pointer :: fieldpg(:,:,:)
    complex(C_FLOAT_COMPLEX), pointer :: fieldpc(:,:,:)

    real sigmachi
    ! DECLARATIONS END
    !=====================================================================

    !=====================================================================
    ! EXCUTABLE BEGIN
  
    ! --------------------------------------------------------------------
    ! GENERATE WHITE NOISE, CONVOLVE WITH P(K), AND OUTPUT
    ! --------------------------------------------------------------------

    delta  => fieldp
    deltac => fieldpc
    deltag => fieldpg 
    
    delta = deltag

    if(fieldid>0.or.fieldid==-1) then

       if(NonGauss == 0) then
          code='rho'
          call generate_noise
          call convolve_noise

       elseif(NonGauss > 0) then
          !1 for classic fNL
          !2 for classic fNL with uncorrelated seed
          !3 for intermittent NG - Gaussian Spike
          !4 for intermittent NG - Chaotic Billiards

          if(NonGauss==1) then 
             call generate_noise             
             code='zetag'
             call convolve_noise
             !now have zeta_lin field. 
             !zeta_NL = zeta_lin + f_NL * [zeta_lin**2 - <zeta_lin**2>]
             lavg = sum(delta(1:n,:,:)**2) / real(n)**3
             call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
             delta = delta + fNL * (delta**2 - avg)
             code='zeta2delta'
             call convolve_noise
             !delta_NL. change code incase fields needs to be written out
             code='rho_fNL'
          
          elseif(NonGauss==2) then
             call generate_noise             
             code='zetag'
             call convolve_noise
             !now have zeta_lin field. 
             !Output field and create uncorrelated zeta**2
             call RandomField_Output
             delta = 0.0
             seed = seed*5
             call generate_noise
             call convolve_noise

             !zeta_NL = zeta_lin + f_NL * [zeta_lin**2 - <zeta_lin**2>]
             delta = delta**2
             lavg = sum(delta(1:n,:,:)) / real(n)**3
             call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
             !get zeta**2 contribution
             delta = fNL * (delta - avg)
             
!             code='zetang'
!             call RandomField_Output
!             code='zetag'
             !add uncorrelated zeta field in place
             call RandomField_Input
             code='zeta2delta'
             call convolve_noise
             !delta_NL. change code incase fields needs to be written out
             code='rho_fNL'

          elseif(NonGauss > 2) then
             code='rho'
             call generate_noise             
             call convolve_noise
             !now have delta_g field. 
             !Output field and create uncorrelated intermittent NG field
             call RandomField_Output
             seed  = seed*5
             delta = 0.0
             code='chi'
             call generate_noise
             call convolve_noise

             if(NonGauss==3) then
                sigmachi = 5.0e-7
                delta = delta/sigmachi !change units of field to sigmas  
             elseif(NonGauss==4) then
                delta = delta + 5.8e-8
             endif
             !now have chi field
             !run through chi to zeta mapping function
             call chi2zetaFNL
             code="zetang"
             call RandomField_Output

             code='zeta2delta'
             call convolve_noise

             code='rho'
             call RandomField_Input
             !delta_NL. change code incase fields needs to be written out
             code='rho_fNL'
          endif
       endif
    endif

    if(fieldid>0.or.fieldid==-2) code='xdisp'
    if(fieldid>0.or.fieldid==-3) code='ydisp'
    if(fieldid>0.or.fieldid==-4) code='zdisp'
    
    if(fieldid<-1) call convolve_noise

    !testing single spherical overdensity
    if(fieldid==-11) then
       !set parameters for profile
       ro = 35.
       rmax = 100. 
       alpha = 1.5
       fcrit = 1.686
       dr = boxsize/n
       ff = 4*3.14159/3 * (rmax/(dr*n))**3
       do i=1,n
          rx = (i-n/2) * dr
          do j=1,n
             ry = (j-n/2) * dr
             do k=1,local_nz
                rz = (k+local_z_start - n/2)*dr
                ar=sqrt(rx**2+ry**2+rz**2+dr)
                if(ar < rmax) then
                   delta(i,j,k) = (3-alpha)/3 * fcrit * (ar/ro)**(-alpha)
                else
                   delta(i,j,k) = ff/(ff-1) * fcrit * (rmax/ro)**(-alpha) 
                endif
             enddo
          enddo
       enddo

       lavg = sum(delta(1:n,:,:)) / real(n)**3
       call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
            mpi_comm_world,ierr)
       lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
       call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
            mpi_comm_world,ierr)
       
       if(myid==0) write(*,*) "avg, sigma of spherical profile = ", avg,sqrt(sig)
          
       call mpi_barrier(mpi_comm_world,ierr)

    endif

    return

  end subroutine RandomField_make

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine RandomField_output

    use intreal_types
    use cosmology
    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use timing_diagnostics
    use globalvars

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    integer ifile
    character *128 outputfile

    real a(1)

    integer iorig
    
    real r4arg
    integer i4arg, indlnb, iargc

    if(myid==0) call timer_begin

    if(code=='rho') then
       outputfile='Fvec_'//outcode
    elseif(code=='rho_fNL') then
       outputfile='Fvec_fNL_'//outcode
    elseif(code=='zetang') then
       outputfile='zetang_'//outcode
    elseif(code=='zetag') then
       outputfile='zetag_'//outcode
    elseif(code=='xdisp') then
       outputfile='etax_'//outcode
    elseif(code=='ydisp') then
       outputfile='etay_'//outcode
    elseif(code=='zdisp') then
       outputfile='etaz_'//outcode
    else
       write(*,*) 'error in code'
       call mpi_finalize(ierr)
       stop
    endif

    open(unit=33,file=outputfile,access='stream')

    offset_tile = n*n*local_z_start
    offset_bytes = offset_tile*int(4,8)+1
    write(33,pos=offset_bytes) (((delta(i,j,k),i=1,n),j=1,n),k=1,local_nz)
    close(33)

    if(myid==0) then
       timing_diagnostics_code = 'writing data'
       call timer_end
    endif
    
    return
    
  end subroutine RandomField_output

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine RandomField_input

    use intreal_types
    use cosmology
    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use timing_diagnostics
    use globalvars

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    integer ifile
    character *128 outputfile

    real a(1)

    integer iorig
    
    real r4arg
    integer i4arg, indlnb, iargc
    real deltaini
    if(myid==0) call timer_begin

    if(code=='rho') then
       outputfile='Fvec_'//outcode
    elseif(code=='rho_fNL') then
       outputfile='Fvec_fNL_'//outcode
    elseif(code=='zetang') then
       outputfile='zetang_'//outcode
    elseif(code=='zetag') then
       outputfile='zetag_'//outcode
    elseif(code=='xdisp') then
       outputfile='etax_'//outcode
    elseif(code=='ydisp') then
       outputfile='etay_'//outcode
    elseif(code=='zdisp') then
       outputfile='etaz_'//outcode
    else
       write(*,*) 'error in code'
       call mpi_finalize(ierr)
       stop
    endif

    open(unit=33,file=outputfile,access='stream')

    offset_tile = n*n*local_z_start
    offset_bytes = offset_tile*int(4,8)+1
    offset_i = 0

    do k=1,local_nz
       do j=1,n
          do i=1,n
             read(33,pos=offset_bytes+offset_i) deltaini             
             delta(i,j,k) = delta(i,j,k)+deltaini
             offset_i = offset_i+4 !4 bytes
          enddo
       enddo
    enddo
    close(33)


    if(myid==0) then
       timing_diagnostics_code = 'writing data'
       call timer_end
    endif
    
    return
    
  end subroutine RandomField_input
  
end module RandomFieldWrapper


