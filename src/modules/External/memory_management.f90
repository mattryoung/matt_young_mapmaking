module memory_management

  use iso_c_binding
  interface

     subroutine ReportMemory(tag) bind(c,name="ReportMemory")

       use iso_c_binding

       character(C_CHAR), intent(in)  :: tag

     end subroutine ReportMemory

     subroutine GetOverhead() bind(c,name="GetOverhead")

     end subroutine GetOverhead

  end interface

  logical report_memory
  
end module memory_management
