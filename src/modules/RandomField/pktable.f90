module pktable

  use globalvars
  implicit none
  character *128 filepktab
  real pktabmin,pktabmax,pktabminl,pktabmaxl,dpktab,dtab
  real, allocatable:: tsav(:), tsav_in(:)
  real, allocatable:: pksav(:), pksav_in(:)
  real, allocatable:: pkchisav(:), pkchisav_in(:)
  integer npktab,itab

contains

  subroutine read_pktable 

    implicit none

    integer i
    real dummy_real

    open(unit=1,file=filepktab)
    read(1,*,end=20) pktabmin
    npktab=1

10  continue

    read(1,*,end=20) pktabmax
    npktab=npktab+1

    goto 10
20  continue

    pktabminl=log10(pktabmin)
    pktabmaxl=log10(pktabmax)
    dpktab=(pktabmaxl-pktabminl)/(npktab-1)

    close(1)       
    if(.not. allocated(tsav))     allocate(tsav(npktab))
    if(.not. allocated(tsav_in))  allocate(tsav_in(npktab)) 
    if(.not. allocated(pksav))    allocate(pksav(npktab))
    if(.not. allocated(pksav_in)) allocate(pksav_in(npktab))
    if(.not. allocated(pkchisav))    allocate(pkchisav(npktab))
    if(.not. allocated(pkchisav_in)) allocate(pkchisav_in(npktab))

    open(unit=1,file=filepktab)

    if(NonGauss==0) then
       do i=1,npktab
          read(1,*) dummy_real,pksav_in(i)
       enddo
       tsav_in     = 0.0
       pkchisav_in = 0.0
    elseif(NonGauss>0) then
       do i=1,npktab
          read(1,*) dummy_real,pksav_in(i),tsav_in(i),pkchisav_in(i)
       enddo
    endif

    close(1)
    return

  end subroutine read_pktable

end module pktable
