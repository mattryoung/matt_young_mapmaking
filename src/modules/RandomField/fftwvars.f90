module fftwvars

  use types
  use mpivars

  !-----------------------------------------------------------------------
  ! FFTW VARIABLES AND CONSTANTS                                  
  !-----------------------------------------------------------------------

  integer(i8b) plan, iplan
  integer local_z_start,local_nz,local_y_start,local_ny,total_local_sizes
  integer kl
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer,parameter:: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  integer,parameter:: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
  integer,parameter:: FFTW_ESTIMATE=0, FFTW_MEASURE=1  
  integer,parameter:: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
  integer,parameter:: FFTW_USE_WISDOM=16
  integer,parameter:: FFTW_THREADSAFE=128
  integer,parameter:: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
  integer,parameter:: FFMPI_COMM_WORLDTW_SCRAMBLED_INPUT=8192
  integer,parameter:: FFTW_SCRAMBLED_OUTPUT=16384

end module fftwvars
