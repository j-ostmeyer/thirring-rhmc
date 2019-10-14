#include "test_utils.fh"
program test_qmrherm_2
  use dwf3d_lib
  use trial
  use vector, only: X
  use qmrherm_module, only: qmrherm, phi0, qmrhprint => printall
  use dirac
  use gammamatrices
  use gforce
  use params
  use comms
  use comms4
  use comms5
  use comms6
  use test_utils
  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  ! common blocks to function

  ! initialise function parameters
  complex(dp) Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp), allocatable :: Phi0_ref(:, :, :, :, :, :)
  complex(dp), allocatable :: Phi0_orig(:, :, :, :, :, :)
  complex(dp) :: x_ref(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: delta_x(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp), allocatable :: delta_Phi0(:, :, :, :, :, :)
  complex(dp) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

  complex(dp) :: sum_delta_x, sum_delta_Phi0
  real(dp) :: max_delta_x, max_delta_Phi0

  integer :: imass, iflag
  real(dp) :: anum(0:ndiag), aden(ndiag)
  real :: res, am
  integer :: itercg

  integer :: i, j, l, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0

#ifdef MPI
  integer, dimension(12) :: reqs_R, reqs_U, reqs_Phi, reqs_Phi0
  integer :: ierr
  call init_MPI
#endif
  qmrhprint = .false.

  allocate (Phi0_ref(kthird, ksizex_l, ksizey_l, ksizet_l, 4, 25))
  allocate (Phi0_orig(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, 25))
  allocate (delta_Phi0(kthird, ksizex_l, ksizey_l, ksizet_l, 4, 25))

  res = 0.1
  am = 0.05
  imass = 3
  iflag = 2

  anum(0) = 0.5
  do i = 1, ndiag
    anum(i) = real(exp(iunit*i*tau/ndiag))
    aden(i) = real(exp(-iunit*0.5*i*tau/ndiag))
  enddo
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
            do l = 1, 25
              Phi0_orig(ithird, ix, iy, it, j, l) = &
                & 1.7*exp(1.0)*exp(iunit*idx*tau/idxmax) + l
            end do
          end do
        end do
      end do
    end do
  end do
#ifdef MPI
  call start_halo_update_5(4, R, 0, reqs_R)
  call start_halo_update_5(4, Phi, 1, reqs_Phi)
  call start_halo_update_6(4, 25, Phi0_orig, 2, reqs_Phi0)
#endif
  do j = 1, 3
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          idx = ip_x*ksizex_l + ix - 1 &
            & + (ip_y*ksizey_l + iy - 1)*ksize &
            & + (ip_t*ksizet_l + it - 1)*ksize*ksize &
            & + (j - 1)*ksize*ksize*ksizet
          u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
          dSdpi(ix, iy, it, j) = real(tau*exp(iunit*idx*tau/idxmax), sp)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(3, u, 3, reqs_u)
  call complete_halo_update(reqs_R)
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_Phi0)
  call complete_halo_update(reqs_u)
#else
  call update_halo_6(4, 25, Phi0_orig)
  call update_halo_5(4, Phi)
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
    Phi0 = Phi0_orig
    call qmrherm(Phi, X, res, itercg, am, imass, anum, aden, ndiag, iflag)
  end do
  ! check output
  if (generate) then
#ifdef MPI
    write_file(x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :), 'test_qmrherm_2_x.dat', MPI_Double_Complex)
    write_file(Phi0(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :), 'test_qmrherm_2_Phi0.dat', MPI_Double_Complex)
#else
    write (6, *) "Generation not possible"
    call exit(1)
#endif
  else
    read_file(x_ref, 'test_qmrherm_2_x.dat', MPI_Double_Complex)
    read_file(Phi0_ref, 'test_qmrherm_2_Phi0.dat', MPI_Double_Complex)

    delta_x = x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - x_ref
    delta_Phi0 = Phi0(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :) - Phi0_ref

    check_equality(itercg, 2, 'itercg', 'test_qmrherm_2')
    check_sum(delta_x, 1e-12, 'x', sum_delta_x, MPI_Double_Complex, 'test_qmrherm_2')
    check_max(delta_x, 1e-14, 'x', max_delta_x, MPI_Double_Precision, 'test_qmrherm_2')
    check_sum(delta_Phi0, 1e-12, 'Phi0', sum_delta_Phi0, MPI_Double_Complex, 'test_qmrherm_2')
    check_max(delta_Phi0, 1e-14, 'Phi0', max_delta_Phi0, MPI_Double_Precision, 'test_qmrherm_2')
  end if
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program test_qmrherm_2
