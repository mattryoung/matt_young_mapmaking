module arrays

  use intreal_types

  integer, parameter :: n1         =  272
  integer, parameter :: n2         =  n1
  integer, parameter :: n3         =  n1
  integer, parameter :: np1d       =  n1
  integer, parameter :: np2d       =  n1
  integer, parameter :: np3d       =  n1
  integer, parameter :: nFmax      =  n1*n2*n3
  integer, parameter :: ircmax     =  max(n1,n2,n3)/2

  integer, parameter :: itabmax    =  2001, &
                        np         =  (np1d)**3, &
	 	        npt        =  (np1d+1)**3
  integer, parameter :: Npv_evmax =   10,   &
                        Npv_evmin  = -10,   &
                        Nevmax     =  19,   &
                        Nfscmax    =  10
  integer, parameter :: nbormax    =  500000
  integer, parameter :: nodemax    =  100000000
  integer, parameter :: ntabmax    =  2000

  integer, parameter :: npkmax=10
  integer, parameter :: npkmaxl=500000
  integer, parameter :: nclmax=30
  integer, parameter :: nboxmax=100000
  integer, parameter :: nhalomax=100

  ! LOCAL GRID VARIABLES
  real(C_FLOAT),            pointer :: delta(:,:,:),delta_u(:,:,:)
  real(C_FLOAT),            pointer :: etax(:,:,:), etay(:,:,:), etaz(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: deltac(:,:,:),deltac_u(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: etaxc(:,:,:), etayc(:,:,:), etazc(:,:,:)

  real(C_FLOAT), pointer :: F(:,:,:)
  type(C_PTR)                       :: deltap, delta_up,etaxp, etayp, etazp
 
  
  
  ! GLOBAL GRID VARIABLES
  real(C_FLOAT),            pointer :: deltag(:,:,:)
  real(C_FLOAT),            pointer :: etaxg(:,:,:), etayg(:,:,:), etazg(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: deltagc(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: etaxgc(:,:,:), etaygc(:,:,:), etazgc(:,:,:)

  type(C_PTR)                       :: deltagp, etaxgp, etaygp, etazgp

  ! CLOUD
  integer, allocatable :: irs2(:)
  integer(kind=2), allocatable :: ixsvec(:,:)
  integer(kind=1), allocatable :: mask(:,:,:)

  ! HALOS
  integer, allocatable :: iwas(:,:),lagrange(:,:)
  integer, allocatable :: iwasa(:,:),lagrangea(:,:)

  real, allocatable :: xpk(:,:),ypk(:,:),zpk(:,:)
  real, allocatable :: vx2pk(:,:),vy2pk(:,:),vz2pk(:,:)
  real, allocatable :: vxpk(:,:),vypk(:,:),vzpk(:,:)
  real, allocatable :: strain_bar(:,:,:)
  real, allocatable :: Fcollv(:,:),RTHLv(:,:),vTHvir(:,:)
  real, allocatable :: hp(:,:),F_ev(:,:),F_pv(:,:)

  real, allocatable :: xpka(:,:),ypka(:,:),zpka(:,:)
  real, allocatable :: vx2pka(:,:),vy2pka(:,:),vz2pka(:,:)
  real, allocatable :: vxpka(:,:),vypka(:,:),vzpka(:,:)
  real, allocatable :: strain_bara(:,:,:)
  real, allocatable :: Fcollva(:,:),RTHLva(:,:),vTHvira(:,:)
  real, allocatable :: hpa(:,:),F_eva(:,:),F_pva(:,:)

  type(C_PTR) :: outbufferp

  ! BOXES
  integer, allocatable :: boxlist(:)
  real, allocatable :: xbox(:),ybox(:),zbox(:)

contains

  subroutine allocate_halos(num_redshifts)

    implicit none 

    integer num_redshifts

    allocate(iwas(npkmaxl,num_redshifts))
    allocate(lagrange(npkmaxl,num_redshifts))
    allocate(xpk(npkmaxl,num_redshifts))
    allocate(ypk(npkmaxl,num_redshifts))
    allocate(zpk(npkmaxl,num_redshifts))
    allocate(vxpk(npkmaxl,num_redshifts))
    allocate(vypk(npkmaxl,num_redshifts))
    allocate(vzpk(npkmaxl,num_redshifts))
    allocate(vx2pk(npkmaxl,num_redshifts))
    allocate(vy2pk(npkmaxl,num_redshifts))
    allocate(vz2pk(npkmaxl,num_redshifts))
    allocate(strain_bar(6,npkmaxl,num_redshifts))
    allocate(Fcollv(npkmaxl,num_redshifts))
    allocate(RTHLv(npkmaxl,num_redshifts))
    allocate(vTHvir(npkmaxl,num_redshifts))
    allocate(hp(npkmaxl,num_redshifts))
    allocate(F_ev(npkmaxl,num_redshifts))
    allocate(F_pv(npkmaxl,num_redshifts))

    allocate(iwasa(npkmax,num_redshifts))
    allocate(lagrangea(npkmax,num_redshifts))
    allocate(xpka(npkmax,num_redshifts))
    allocate(ypka(npkmax,num_redshifts))
    allocate(zpka(npkmax,num_redshifts))
    allocate(vxpka(npkmax,num_redshifts))
    allocate(vypka(npkmax,num_redshifts))
    allocate(vzpka(npkmax,num_redshifts))
    allocate(vx2pka(npkmax,num_redshifts))
    allocate(vy2pka(npkmax,num_redshifts))
    allocate(vz2pka(npkmax,num_redshifts))
    allocate(strain_bara(6,npkmax,num_redshifts))
    allocate(Fcollva(npkmax,num_redshifts))
    allocate(RTHLva(npkmax,num_redshifts))
    allocate(vTHvira(npkmax,num_redshifts))
    allocate(hpa(npkmax,num_redshifts))
    allocate(F_eva(npkmax,num_redshifts))
    allocate(F_pva(npkmax,num_redshifts))

    RTHLv=0

  end subroutine allocate_halos

  subroutine allocate_boxes

    allocate(xbox(nboxmax))
    allocate(ybox(nboxmax))
    allocate(zbox(nboxmax))
    allocate(boxlist(nboxmax))

  end subroutine allocate_boxes

end module arrays
