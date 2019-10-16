module params
  implicit none
  save

  ! Type definitions
  integer, parameter :: k4b = selected_int_kind(9)
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: sp = kind(1.)

  ! Lattice parameters
#define KSIZE 12
#define KSIZET 12
  integer, parameter :: ksize = KSIZE, ksizet = KSIZET
  integer, parameter :: kthird = 24
  integer, parameter :: kvol = ksize*ksize*ksizet
  integer, parameter :: ndiag = 25, ndiagg = 12
  integer, parameter :: Nf = 1
  real(dp), parameter :: akappa = 0.5d0
#ifndef MPI
  integer, parameter :: ksizex_l = ksize, ksizey_l = ksize, ksizet_l = ksizet
  integer, parameter :: kvol_l = kvol
  integer, parameter :: np_x = 1, np_y = 1, np_t = 1, np_global = 1
  integer, parameter :: ip_x = 0, ip_y = 0, ip_t = 0, ip_global = 0
#else
#if !(defined(NP_X) && defined(NP_Y) && defined(NP_T))
#error "NP_X, NP_Y, and NP_T must be defined for MPI compilation."
#endif
#if (KSIZE / NP_X) * NP_X != KSIZE
#error "ksize must be divisible by NP_X"
#elif (KSIZE / NP_Y) * NP_Y != KSIZE
#error "ksize must be divisible by NP_Y"
#elif (KSIZET / NP_T) * NP_T != KSIZET
#error "ksizet must be divisible by NP_T"
#endif
  integer, parameter :: np_x = NP_X, np_y = NP_Y, np_t = NP_T
  integer, parameter :: ksizex_l = ksize/np_x
  integer, parameter :: ksizey_l = ksize/np_y
  integer, parameter :: ksizet_l = ksizet/np_t
#endif

  ! Control parameters
  integer, parameter :: istart = -1
  integer, parameter :: iread = 1
  integer, parameter :: iwrite = 0
  integer, parameter :: iprint = 5
  integer, parameter :: iseed = 1
  integer, parameter :: icheckpoint = 100

  ! Inverter
  integer :: max_qmr_iters = 2   !QMRHERM
  integer :: niterc = kthird*kvol ! CONGRAD

  ! inverter residuals
  real, parameter :: respbp = 1.0e-6, rescgg = 1.0e-6
  real, parameter :: rescga = 1e-9
  real, parameter :: rescgm = 1e-9
  logical, parameter :: spmd = .false.

  ! max step in molecular dynamics evolution
  integer, parameter :: itermax = 1000

  ! Runtime parameters
  real :: beta
  real :: am3
  integer :: ibound
end module params
