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
    
    integer count,count_rate

!    call cpu_time(timer_real)
!    timing_diagnostics_t1=timer_real
    call system_clock(count,count_rate)    
    timing_diagnostics_t1=count
    
  END SUBROUTINE timer_begin
  
  ! =====================================================

  SUBROUTINE timer_end
    
    use textlib 
    
    implicit none
    
    integer count,count_rate
    real dt
    
!    call cpu_time(timer_real)
!    timing_diagnostics_t2=timer_real
    call system_clock(count,count_rate)    
    timing_diagnostics_t2=count
    dt = timing_diagnostics_t2 - timing_diagnostics_t1
    dt = dt / count_rate
    
    write(*,19) &
         timing_diagnostics_code(1:indlnb(timing_diagnostics_code)),dt
19  format(/,1x,a51,': ',f6.1,' sec ')
    
  END SUBROUTINE timer_end

END MODULE timing_diagnostics

