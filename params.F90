module params
  implicit none
  save

  ! Type definitions
  integer, parameter :: k4b=selected_int_kind(9)
  integer, parameter :: dp=kind(1.d0)

  ! Lattice parameters
#define KSIZE 12
#define KSIZET 12
  integer, parameter :: ksize=KSIZE, ksizet=KSIZET
  integer, parameter :: kthird=4
  integer, parameter :: kvol=ksize*ksize*ksizet
  integer, parameter :: ndiag=25, ndiagg=12
  integer, parameter :: Nf=1
  real, parameter :: akappa = 0.5
#ifndef MPI
  integer, parameter :: ksizex_l=ksize, ksizey_l=ksize, ksizet_l=ksizet
  integer, parameter :: kvol_l = kvol
  integer, parameter :: np_x=1, np_y=1, np_t=1, np_global=1
  integer, parameter :: ip_x=0, ip_y=0, ip_t=0, ip_global=0
#else
#if !(defined(NP_X) && defined(NP_Y) && defined(NP_T))
#error "NP_X, NP_Y, and NP_T must be defined for MPI compilation."
#endif
#if (ksize / NP_X) * NP_X != ksize
#error "ksize must be divisible by NP_X"
#elif (ksize / NP_Y) * NP_Y != ksize
#error "ksize must be divisible by NP_Y"
#elif (ksizet / NP_T) * NP_T != ksizet
#error "ksizet must be divisible by NP_T"
#endif
  integer, parameter :: np_x=NP_X, np_y=NP_Y, np_t=NP_T
  integer, parameter :: ksizex_l = ksize / np_x
  integer, parameter :: ksizey_l = ksize / np_y
  integer, parameter :: ksizet_l = ksizet / np_t
#endif
  
  ! Runtime parameters
  real :: beta
  real :: am3
  integer :: ibound
end module params
