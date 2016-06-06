module gaussian_field

  use grid
  use fftw_interface
  use mpivars
  use timing_diagnostics
  use globalvars
  use openmpvars

  implicit none

  !-----------------------------------------------------------------------
  ! BOOKKEEPING AND GLOBAL VARIABLES
  !-----------------------------------------------------------------------

  real xflag,yflag,zflag,rhoflag
  real kx,ky,kz,kmax,ak,dk,d3k,dq
  real z
  real D

  integer seed

contains

  subroutine generate_noise

    use random

    implicit none

    real*8 lavg,avg,lsig,sig
    real work(1)
    real(dp) :: xr 

    integer :: seeds(MaxRandNumStreams,IRandNumSize)   

    if(myid==0.and.rf_report==1) call timer_begin
    ! White noise
    call rans(ntasks,seed,seeds)
    do k=1,local_nz
       do j=1,n
          do i=1,n
             call gaussdev(seeds(myid+1,:),xr)
             delta(i,j,k)=xr
          enddo
       enddo
    enddo

    ! Report mean and variance
    lavg = sum(delta(1:n,:,:)) / real(n)**3
    call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
    lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
    call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)

    if(myid==0.and.rf_report==1) then
       write(*,96) avg,sqrt(sig)
       timing_diagnostics_code='white noise'
       call timer_end
    endif

    return
96  format(3x,'white noise:',       /5x' mean = ',&
           1pe13.6,/5x,'sigma = ',1pe13.6)
    
  end subroutine generate_noise
  
  !=======================================================================
  
  subroutine convolve_noise

    use pktable

    implicit none    

    complex(C_FLOAT_COMPLEX) ctemp

    real, allocatable :: work(:)

    double precision lavg,avg,lsig,sig

    real dummy_real
    real tf
    real pk
    real pkchi
!    pksav=sqrt(pksav_in*d3k*n**3)
!    tsav=sqrt(tsav_in*d3k*n**3)

    ! --------------------------------------------------------------------
    ! FFT FROM REAL SPACE TO K SPACE
    ! --------------------------------------------------------------------

    if(myid==0.and.rf_report==1) call timer_begin

    call fftwf_mpi_execute_dft_r2c(plan, delta, deltac)

    if(myid==0.and.rf_report==1) then
       timing_diagnostics_code='real2complex'
       call timer_end
    endif
    
    ! --------------------------------------------------------------------
    ! MULTIPLY delta(k1,k2,k3) by P(k)                                    
    ! --------------------------------------------------------------------
    
    if(myid==0.and.rf_report==1) call timer_begin

!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(deltac) 
    do k=1,local_nz
       kz=(k+local_z_start-1)*dk
       if (k+local_z_start.gt.n12) kz=kz-kmax
       
       do j=1,n
          ky=(j-1)*dk
          if (j.gt.n12) ky=ky-kmax
          
          do i=1,n12+1
             kx=(i-1)*dk
             if (i.gt.n12) kx=kx-kmax
             
             ak=sqrt(kx**2+ky**2+kz**2)
             
             itab=int((log10(ak)-pktabminl)/dpktab)+1
             if(ak==0.) then
                pk    = 0.
                pkchi = 0.
                tf    = 0.
             else                                                
                dtab=(log10(ak)-(pktabminl+(itab-1)*dpktab))/dpktab
                pk    = pksav(itab)*(1.-dtab)+pksav(itab+1)*dtab
                tf    = tsav(itab)*(1.-dtab)+tsav(itab+1)*dtab
                pkchi = pkchisav(itab)*(1.-dtab)+pkchisav(itab+1)*dtab
             endif
             
             ctemp = deltac(i,j,k)
             
             if(code.eq.'rho') then
                ctemp = ctemp * pk 
             elseif(code.eq.'chi') then
                ctemp = ctemp * pkchi 
             elseif(code.eq.'zetag') then
                ctemp = ctemp * pk/tf**2 
             elseif(code.eq.'zeta2delta') then
                ctemp = ctemp * tf**2 
             elseif(code.eq.'xdisp') then
                ctemp = -ctemp * kx / ak**2 * cmplx(0.0,1.0)
                if(i.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'ydisp') then
                ctemp = -ctemp * ky / ak**2 * cmplx(0.0,1.0)
                if(j.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'zdisp') then
                ctemp = -ctemp * kz / ak**2 * cmplx(0.0,1.0)
                if(k+local_z_start.eq.n12+1) ctemp = 0.0
             endif
             
             if(ak==0.) ctemp = 0.0

             deltac(i,j,k) = ctemp
             
          enddo
       enddo
    enddo
    
    if(myid==0.and.rf_report==1) then
       timing_diagnostics_code='convolution'
       call timer_end
    endif

    ! --------------------------------------------------------------------
    ! FFT FROM K SPACE TO REAL SPACE
    ! --------------------------------------------------------------------
    
    if(myid==0.and.rf_report==1) call timer_begin
    
    call fftwf_mpi_execute_dft_c2r(iplan,deltac,delta)

    if(myid==0.and.rf_report==1) then
       timing_diagnostics_code='complex2real'
       call timer_end
    endif
   
    ! Normalization of backwards FFT
    do k=1,local_nz
        do j=1,n
            do i=1,n
                delta(i,j,k) = delta(i,j,k) / real(n)**3
            enddo
        enddo
    enddo
    ! Report mean and variance
    lavg = sum(delta(1:n,:,:)) / real(n)**3
    call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
    lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
    call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
    
    if(myid==0.and.rf_report==1) then
       if(code.eq.'rho')         write(*,96) avg,sqrt(sig)
       if(code.eq.'xdisp')       write(*,97) avg,sqrt(sig)
       if(code.eq.'ydisp')       write(*,98) avg,sqrt(sig)
       if(code.eq.'zdisp')       write(*,99) avg,sqrt(sig)
       if(code.eq.'zetag')       write(*,100) avg,sqrt(sig)
       if(code.eq.'zeta2delta')  write(*,101) avg,sqrt(sig)
       if(code.eq.'chi')         write(*,102) avg,sqrt(sig)
    endif


96  format(3x,'density:',       /5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
97  format(3x,'x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
98  format(3x,'y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
99  format(3x,'z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
100  format(3x,'zeta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
101  format(3x,'zeta2delta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
102  format(3x,'chi:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)

    call mpi_barrier(mpi_comm_world,ierr)
    
    return
    
  end subroutine convolve_noise

end module gaussian_field
