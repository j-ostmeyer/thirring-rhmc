#include "test_utils.fh"
program test_congrad
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
  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real*8, parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  complex(dp) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: x_ref(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: diff(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: sum_diff
  real(dp) :: max_diff

  integer :: imass, iflag, isweep, iter
  real :: res, am
  integer :: itercg

  integer :: i, j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0
#ifdef MPI
  integer, dimension(12) :: reqs_X, reqs_Phi, reqs_u
  integer :: ierr

  call init_MPI
  !call gdbwait
#endif

  res = 0.1
  am = 0.05
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
  do i = 1, timing_loops
    call congrad(Phi, res, itercg, am, imass)
  end do
  ! check output
  if (generate) then
    write_file(x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :), 'test_congrad.dat', MPI_Double_Complex)
  else
    read_file(x_ref, 'test_congrad.dat', MPI_Double_Complex)

    diff = x_ref - x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
    sum_diff = sum(diff)
    max_diff = maxval(abs(diff))
#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, sum_diff, 1, MPI_Double_Complex, MPI_Sum, comm, ierr)
    call MPI_AllReduce(MPI_IN_PLACE, max_diff, 1, MPI_Double_Precision, MPI_Max, comm, ierr)
#endif
    if (ip_global .eq. 0) then
      if (itercg .ne. 27) then
        print *, 'itercg looks wrong: ', itercg, ' != 27'
      end if
      if (abs(sum_diff) .gt. 2) then
        print *, 'sum delta too large: ', sum_diff
      end if
      if (max_diff .gt. 5e-2) then
        print *, 'max delta too large: ', max_diff
      end if
    end if
  end if
  call MPI_Finalize(ierr)
end program
