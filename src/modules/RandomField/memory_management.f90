module memory_management

  use fftwvars
  use mpivars  

  implicit none

  character *12 memtag
  real allocated
  integer gridsize

contains

  subroutine report_allocation(allocatecur)

    implicit none

    real allocatecur

    allocatecur = allocatecur / 1.e9
    allocated = allocated + allocatecur

    if(myid==0) then
       write(*,*) '*************** ALLOCATE ******************'
       write(*,*) allocatecur,' GB allocated for ',memtag
       write(*,*) '-------------------------------------------'
       write(*,*) ' TOTAL MEMORY ALLOCATED: ',allocated, 'GB'
       write(*,*) '       ',allocated/gridsize**2/local_nz*1e9,' bytes/cell'
       write(*,*) '*******************************************'
    endif
    
    return
    
  end subroutine report_allocation

  subroutine report_deallocation(allocatecur)

    implicit none

    real allocatecur

    allocatecur = allocatecur / 1.e9 
    allocated = allocated - allocatecur 

    if(myid==0) then
       write(*,*) '************** DEALLOCATE *****************'
       write(*,*) allocatecur,' GB deallocated for ',memtag
       write(*,*) '-------------------------------------------'
       write(*,*) ' TOTAL MEMORY ALLOCATED: ',allocated, 'GB'
       write(*,*) '       ',allocated/gridsize**2/local_nz*1e9,' bytes/cell'
       write(*,*) '*******************************************'
    endif

    return

  end subroutine report_deallocation

end module memory_management
