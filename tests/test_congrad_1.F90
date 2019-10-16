#include "test_utils.fh"
program test_congrad_1
  use dwf3d_lib
  use trial, only: u
  use vector
  use dirac
  use comms
  use comms4
  use comms5
  use gammamatrices
  use measure_module
  use test_utils
  use inverter_utils
  implicit none

  ! general parameters
  complex, parameter :: iunit = cmplx(0, 1)
  real*8, parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  complex(dp) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp), allocatable :: delta_Phi(:, :, :, :, :)
  real(dp) :: sum_delta_Phi

  integer :: imass, iflag, isweep, iter
  real :: res, am
  real(dp) :: h, hg, hp, s
  integer :: itercg

  integer :: j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0
#ifdef MPI
  integer, dimension(12) :: reqs_X, reqs_Phi, reqs_u
  integer :: ierr
  call init_MPI
#endif
  allocate (delta_Phi(kthird, ksizex_l, ksizey_l, ksizet_l, 4))

  h = 0
  hg = 0
  hp = 0
  s = 0
  res = 1.0e-6
  am = 0.2
  imass = 3
  iflag = 0
  isweep = 1
  iter = 0

  do j = 1, 4
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          do ithird = 1, kthird
            idx = ithird + (ip_x*ksizex_l + ix - 1)*kthird &
              & + (ip_y*ksizey_l + iy - 1)*kthird*ksize &
              & + (ip_t*ksizet_l + it - 1)*kthird*ksize*ksize &
              & + (j - 1)*kthird*ksize*ksize*ksizet
            Phi(ithird, ix, iy, it, j) = 1.1*exp(iunit*idx*tau/idxmax)
            X(ithird, ix, iy, it, j) = 0.5*exp(1.0)*exp(iunit*idx*tau/idxmax)
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, X, 0, reqs_X)
  call start_halo_update_5(4, Phi, 0, reqs_Phi)
#endif
  idx = 0
  do j = 1, 3
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          idx = ip_x*ksizex_l + ix &
            & + (ip_y*ksizey_l + iy - 1)*ksize &
            & + (ip_t*ksizet_l + it - 1)*ksize*ksize &
            & + (j - 1)*ksize*ksize*ksizet
          u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(3, u, 1, reqs_u)
  call complete_halo_update(reqs_X)
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_u)
#else
  call update_halo_5(4, Phi)
  call update_halo_5(4, X)
  call update_halo_4(3, u)
#endif
  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1

  call init_gammas()
  ! call function
  call congrad(Phi, res, itercg, am, imass) ! results are in vector.x
  ! check convergence
  call dirac_operator(xout, x, u, am, imass)
  delta_Phi = Phi(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - &
    &                xout(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
  delta_Phi = abs(delta_Phi)**2

  check_sum(delta_Phi, 1e-6, 'xout', sum_delta_Phi, MPI_Double_Precision, 'test_congrad_1')

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program

