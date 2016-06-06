program ffttest

  use, intrinsic :: iso_c_binding

  implicit none

  include 'mpif.h'

  include 'fftw3-mpi.f03'
  integer(C_INT) i,j,k,n

  integer ierr, myid, ntasks, request
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status
  integer, allocatable :: requests(:)

  type(C_PTR) :: plan, iplan, p
  integer(C_INTPTR_T) :: n0,total_local_sizes, local_nz, local_z_start
  real(C_FLOAT), pointer :: data(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: datac(:,:,:)
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  n = 480
  n0 = n
  total_local_sizes = fftwf_mpi_local_size_3d(n0,n0,n0,mpi_comm_world,&
                                             local_nz,local_z_start)

  p = fftwf_alloc_complex(total_local_sizes)  
  call c_f_pointer(p, data, [2*(n0/2+1),n0,local_nz])
  call c_f_pointer(p,datac, [  (n0/2+1),n0,local_nz])

  plan  = fftwf_mpi_plan_dft_r2c_3d(n0,n0,n0, data,datac,mpi_comm_world,&
                                    fftw_estimate)
  iplan = fftwf_mpi_plan_dft_c2r_3d(n0,n0,n0,datac, data,mpi_comm_world,&
                                    fftw_estimate)
  do k=1,local_nz
     do j=1,n
        do i=1,n
           data(i,j,k)=i+j+k
        enddo
     enddo
  enddo

  call fftwf_mpi_execute_dft_r2c( plan, data, datac)
  call fftwf_mpi_execute_dft_c2r(iplan,datac, data)

end program ffttest
