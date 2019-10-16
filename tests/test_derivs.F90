#include "test_utils.fh"
program test_derivs
  use dwf3d_lib
  use dirac
  use gforce
  use comms
  use comms4
  use comms5
  use test_utils
  use derivs_module, only: derivs
  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  ! common blocks to function

  ! initialise function parameters
  complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  complex(dp) Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)
  complex(dp) R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real :: diff(ksizex_l, ksizey_l, ksizet_l, 3), sum_diff, max_diff

  real(dp), parameter :: anum = tau
  integer, parameter :: iflag = 0

  integer :: i, j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0
#ifdef MPI
  integer, dimension(12) :: reqs_R, reqs_Phi, reqs_X2
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
            X2(ithird, ix, iy, it, j) = 0.5*exp(1.0)*exp(iunit*idx*tau/idxmax)
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, R, 0, reqs_R)
  call start_halo_update_5(4, Phi, 0, reqs_Phi)
  call start_halo_update_5(4, X2, 0, reqs_X2)
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
          dSdpi_ref(ix, iy, it, j) = real(tau*exp(iunit*idx*tau/idxmax), sp)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call complete_halo_update(reqs_R)
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_X2)
#else
  call update_halo_5(4, X2)
  call update_halo_5(4, R)
  call update_halo_4(3, u)
#endif
  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1

  call init_gammas()
  ! call function
  do i = 1, timing_loops
    dSdpi = dSdpi_ref
    call derivs(R, X2, anum, iflag)
  end do
  ! check output
  if (generate) then
    write_file(dSdpi, 'test_derivs.dat', MPI_Real)
  else
    read_file(dSdpi_ref, 'test_derivs.dat', MPI_Real)

    diff = dSdpi - dSdpi_ref
    sum_diff = sum(diff)
    max_diff = maxval(abs(diff))
#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, sum_diff, 1, MPI_Real, MPI_Sum, &
      & comm, ierr)
    call MPI_AllReduce(MPI_IN_PLACE, max_diff, 1, MPI_Real, MPI_Max, &
      & comm, ierr)
#endif
    if (ip_global .eq. 0) then
      if (abs(sum_diff) .gt. 0.3) then
        print *, 'sum delta too large: ', sum_diff
      end if
      if (max_diff .gt. 0.01) then
        print *, 'max delta too large: ', max_diff
      end if
    end if
  end if
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program test_derivs
