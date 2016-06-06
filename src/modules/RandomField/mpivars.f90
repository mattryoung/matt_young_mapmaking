module mpivars

  !-----------------------------------------------------------------------
  ! MPI AND FFTW VARIABLES AND CONSTANTS                                  
  !-----------------------------------------------------------------------

  include 'mpif.h'

  integer ierr, myid, ntasks
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status

end module mpivars
