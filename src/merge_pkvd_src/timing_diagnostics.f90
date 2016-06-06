MODULE timing_diagnostics

  ! A module to time and report sections of code

  implicit none
  
  double precision timing_diagnostics_t1,timing_diagnostics_t2
  character*128 timing_diagnostics_code 
  real timer_real
  real, allocatable :: section_timers(:)
  
CONTAINS
  
  ! =====================================================

  SUBROUTINE timer_begin
    
    implicit none
    
    call cpu_time(timer_real)
    timing_diagnostics_t1=timer_real
    
  END SUBROUTINE timer_begin
  
  ! =====================================================

  SUBROUTINE timer_end
    
    use textlib 
    
    implicit none
    
    integer nt
    real dt
    
    call cpu_time(timer_real)
    timing_diagnostics_t2=timer_real

    dt = timing_diagnostics_t2 - timing_diagnostics_t1
    
    write(*,19) &
         timing_diagnostics_code(1:indlnb(timing_diagnostics_code)),dt
19  format(/,1x,a51,': ',f6.1,' sec ')
    
  END SUBROUTINE timer_end

END MODULE timing_diagnostics

