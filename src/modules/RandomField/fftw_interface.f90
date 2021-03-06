module fftw_interface

  use, intrinsic :: iso_c_binding
  use mpivars
  use openmpvars

  implicit none

  ! MPI
  include 'fftw3-mpi.f03'
  type(C_PTR)         :: plan, iplan
  integer(C_INTPTR_T) :: nfft,total_local_sizes, local_nz,local_z_start

  ! SERIAL
  type(C_PTR)         :: plan_s, iplan_s
  integer(C_INT) :: nfft_s
 
  ! ARRAY POINTERS FOR PLAN MAKING
  real(C_FLOAT),            pointer :: fftw_in(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: fftw_out(:,:,:)

contains

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

  subroutine fftw_allocate_parallel(ext_array,ext_arrayc,ext_pointer)
    
    implicit none

    type(C_PTR)                       :: ext_pointer
    real(C_FLOAT),            pointer :: ext_array(:,:,:)
    complex(C_FLOAT_COMPLEX), pointer :: ext_arrayc(:,:,:)
    integer(C_SIZE_T)                    fftw_alloc_size

    fftw_alloc_size = 2*(nfft/2+1)*nfft*local_nz

    ext_pointer   = fftwf_alloc_real(fftw_alloc_size)
    call c_f_pointer(ext_pointer, ext_array, [2*(nfft/2+1),nfft,local_nz])
    call c_f_pointer(ext_pointer,ext_arrayc, [  (nfft/2+1),nfft,local_nz])

    return 

  end subroutine fftw_allocate_parallel
    
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

  subroutine fftw_allocate_serial(ext_array,ext_arrayc,ext_pointer)
    
    implicit none

    type(C_PTR)                       :: ext_pointer
    real(C_FLOAT),            pointer :: ext_array(:,:,:)
    complex(C_FLOAT_COMPLEX), pointer :: ext_arrayc(:,:,:)
    integer(C_SIZE_T)                    fftw_alloc_size

    fftw_alloc_size = (nfft_s+2)*nfft_s**2
    ext_pointer   = fftwf_alloc_real(fftw_alloc_size)
    call c_f_pointer(ext_pointer, ext_array, [2*(nfft_s/2+1),nfft_s,nfft_s])
    call c_f_pointer(ext_pointer,ext_arrayc, [  (nfft_s/2+1),nfft_s,nfft_s])

  end subroutine fftw_allocate_serial

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

  subroutine fftw_initialize()                            

    implicit none

    ierr = fftwf_init_threads()
    call fftwf_plan_with_nthreads(omp_get_max_threads())

    total_local_sizes = fftwf_mpi_local_size_3d(nfft,nfft,nfft,mpi_comm_world,&
                                             local_nz,local_z_start)

    plan  = fftwf_mpi_plan_dft_r2c_3d(nfft,nfft,nfft,fftw_in ,fftw_out,&
                                      mpi_comm_world,fftw_estimate)
    iplan = fftwf_mpi_plan_dft_c2r_3d(nfft,nfft,nfft,fftw_out,fftw_in ,&
                                      mpi_comm_world,fftw_estimate)

    plan_s  = fftwf_plan_dft_r2c_3d(nfft_s,nfft_s,nfft_s,fftw_in ,fftw_out,&
                                      fftw_estimate)
    iplan_s = fftwf_plan_dft_c2r_3d(nfft_s,nfft_s,nfft_s,fftw_out,fftw_in ,&
                                      fftw_estimate)

    if(myid==0) write(*,101) nfft,local_nz
101 format(/,3x,'Slab decomposition: n = ',i0,' local_nz = ',i0)

    return

  end subroutine fftw_initialize

end module fftw_interface


