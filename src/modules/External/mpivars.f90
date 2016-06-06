module mpivars

  !-----------------------------------------------------------------------
  ! MPI AND FFTW VARIABLES AND CONSTANTS                                  
  !-----------------------------------------------------------------------

  include 'mpif.h'

  integer ierr, myid, ntasks, request
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status
  integer, allocatable :: requests(:)

end module mpivars
