program RandomFieldDriver

  !=======================================================================
  !
  ! GENERATES UNIGRID INITIAL CONDITIONS FOR COSMOLOGICAL SIMULATIONS
  !
  !                                                AUTHOR: MARCELO ALVAREZ
  !                                                LAST EDIT:     10.22.14
  !
  !=======================================================================

  !=======================================================================
  ! DECLARATIONS BEGIN  

  !-----------------------------------------------------------------------
  ! INCLUDE NECESSARY MODULES
  !-----------------------------------------------------------------------

  use mpivars    
  use textlib
  use RandomFieldWrapper

  !-----------------------------------------------------------------------
  ! IMPLICIT STATEMENT
  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------
  ! COMMANDLINE PARAMETERS
  !-----------------------------------------------------------------------

  character *128 filepktab, outcode
  real           boxsize, ainit
  integer        nmesh, nbuff, ntile, seed, fieldid

  ! DATE AND TIME
  character(8)          :: dnt_date
  character(10)         :: dnt_time
  character(5)          :: dnt_zone
  integer, dimension(8) :: dnt_values

  ! DECLARATIONS END
  !=======================================================================

  !=======================================================================
  ! EXCUTABLE BEGIN
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  !-----------------------------------------------------------------------
  ! GET COMMANDLINE PARAMETERS
  ! ----------------------------------------------------------------------

  if(iargc().lt.1) then
     if(myid==0) write(*,99)
99   format('usage: RandomFieldDriver <filepktab> <outcode> [<Lbox [Mpc]> <nmesh> <nbuff> <ntile> <ainit> <seed> <fieldflag>]')
     call mpi_finalize(ierr)  
     stop 
  endif
     
  call getarg(1,filepktab)
  call getarg(2,outcode)
  boxsize = r4arg(3,1.e2)
  nmesh   = i4arg(4,256)
  nbuff   = i4arg(5,16)
  ntile   = i4arg(6,10)
  ainit   = r4arg(7,1e0)
  seed    = i4arg(8,13857)
  fieldid = i4arg(9,1)

  call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
  if(myid==0) write(*,11) dnt_values(3),dnt_values(2),&
       dnt_values(1),dnt_values(5),dnt_values(6),&
       dnt_values(7),dnt_zone

  call RandomField_init(filepktab,outcode,boxsize,nmesh,nbuff,ntile,ainit,seed)       
  call RandomField_make(-1)
  call RandomField_output
  call RandomField_make(-2)
  call RandomField_output
  call RandomField_make(-3)
  call RandomField_output
  call RandomField_make(-4)
  call RandomField_output

  call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
  if(myid==0) write(*,12) dnt_values(3),dnt_values(2),&
       dnt_values(1),dnt_values(5),dnt_values(6),&
       dnt_values(7),dnt_zone

  call mpi_finalize(ierr)
  stop

 11 format(/,3x,61('-'),/,3x,'Random field running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)

 12 format(/,3x,61('-'),/,3x,&
         'Random field finished running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)

end program RandomFieldDriver

