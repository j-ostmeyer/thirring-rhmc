#include "test_utils.fh"
program test_dslashd
  use params
  use dwf3d_lib
  use dirac
  use comms
  use comms4
  use comms5
  use test_utils
  implicit none

  ! general parameters
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  ! common blocks to function

  ! initialise function parameters
  complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  complex(dp) Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) Phiref(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) diff(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) sum_diff
  real(dp) max_diff

  real, parameter :: am = 0.05
  integer, parameter :: imass = 5

  integer :: j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx
#ifdef MPI
  integer, dimension(12) :: reqs_R, reqs_U, reqs_Phi
  integer :: ierr
  call init_MPI
#endif
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
            R(ithird, ix, iy, it, j) = 1.3*exp(iunit*idx*tau/idxmax)
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, R, 0, reqs_R)
  call start_halo_update_5(4, Phi, 1, reqs_Phi)
#endif
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
  !call complete_halo_update(reqs_R) ! done in dslashd
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_u)
#else
  call update_halo_5(4, R)
  call update_halo_5(4, Phi)
  call update_halo_4(3, u)
#endif

  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1

  call init_gammas()
  ! call function
#ifdef MPI
  call dslashd(Phi, R, u, am, imass, reqs_R)
  call start_halo_update_5(4, Phi, 2, reqs_Phi)
  call complete_halo_update(reqs_Phi)
#else
  call dslashd(Phi, R, u, am, imass)
  call update_halo_5(4, Phi)
#endif
  ! check output
  read_file(Phiref, 'test_dslashd_5.dat', MPI_Double_Complex)

  diff = Phi(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - Phiref
  check_max(diff, 1e-11, 'Phi', max_diff, MPI_Double_Precision, 'test_dslashd_5')
  check_sum(diff, 1e-11, 'Phi', sum_diff, MPI_Double_Complex, 'test_dslashd_5')
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
