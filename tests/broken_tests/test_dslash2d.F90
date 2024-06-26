#include "test_utils.fh"
program test_dslash2d
  use params
  use dwf3d_lib
  use dirac
  use comms
  use comms4
  use test_utils
  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  ! common blocks to function

  ! initialise function parameters
  complex(dp) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  complex(dp) :: Phi(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: Phiref(ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: R(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: diff(ksizex_l, ksizey_l, ksizet_l, 4), sum_diff
  real(dp) :: max_diff

  real, parameter :: am = 0.05
  integer, parameter :: imass = 3

  integer :: i, j, ix, iy, it
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet
  integer :: idx
#ifdef MPI
  integer, dimension(12) :: reqs_R, reqs_Phi, reqs_u
  integer :: ierr
  call init_MPI
#endif
  do j = 1, 4
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          idx = ip_x*ksizex_l + ix - 1 &
            & + (ip_y*ksizey_l + iy - 1)*ksize &
            & + (ip_t*ksizet_l + it - 1)*ksize*ksize &
            & + (j - 1)*ksize*ksize*ksizet
          Phi(ix, iy, it, j) = 1.1*exp(iunit*idx*tau/idxmax)
          R(ix, iy, it, j) = 1.3*exp(iunit*idx*tau/idxmax)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(4, Phi, 0, reqs_Phi)
  call start_halo_update_4(4, R, 1, reqs_R)
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
  call start_halo_update_4(3, u, 2, reqs_u)
  call MPI_Waitall(12,reqs_u,MPI_STATUSES_IGNORE,ierr)
  call MPI_Waitall(12,reqs_R,MPI_STATUSES_IGNORE,ierr)
  call MPI_Waitall(12,reqs_Phi,MPI_STATUSES_IGNORE,ierr)
#else
  call update_halo_4(4, Phi)
  call update_halo_4(4, R)
  call update_halo_4(3, u)
#endif

  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1

  call init_gammas()
  ! call function
  do i = 1, timing_loops
    call dslash2d(Phi, R, u)
#ifdef MPI
    call start_halo_update_4(4, Phi, 2, reqs_Phi)
    call MPI_Waitall(12,reqs_Phi,MPI_STATUSES_IGNORE,ierr)
#else
    call update_halo_4(4, Phi)
#endif
  end do
  ! check output
  if (generate) then
    write_file(Phi(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :), 'test_dslash2d.dat', MPI_Double_Complex)
  else
    read_file(Phiref, 'test_dslash2d.dat', MPI_Double_Complex)

    diff = Phi(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - Phiref
    sum_diff = sum(diff)
    max_diff = maxval(abs(diff))
#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, sum_diff, 1, MPI_Double_Complex, MPI_Sum, &
      & comm, ierr)
    call MPI_AllReduce(MPI_IN_PLACE, max_diff, 1, MPI_Double_Precision, MPI_Max, &
      & comm, ierr)
#endif
    check_max(diff, 1e-13, 'Phi', max_diff, MPI_Double_Precision, 'test_dslash2d')
    check_sum(diff, 1e-11, 'Phi', sum_diff, MPI_Double_Complex, 'test_dslash2d')
  end if
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
