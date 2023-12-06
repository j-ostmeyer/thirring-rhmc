#include "test_utils.fh"
program test_qmrherm_0
  ! use dwf3d_lib
  use trial, only: u
  ! use gdbhook
  use qmrherm_module, only: qmrherm, phi0, qmrhprint => printall
  !use dirac
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
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: Phi0_ref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, 25)
  complex(dp) :: Phi0_orig(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, 25)
  complex(dp) :: X2_ref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: delta_x(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: delta_Phi0(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, 25)
  complex(dp) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)

  complex(dp) :: sum_delta_x, sum_delta_Phi0
  real(dp) :: max_delta_x, max_delta_Phi0

  integer :: imass, iflag
  real(dp) :: anum(0:ndiag), aden(ndiag)
  real :: res, am
  integer :: itercg

  integer :: i

#ifdef MPI
  integer, dimension(16) :: reqs_R, reqs_Phi, reqs_X
  integer :: ierr
  call init_MPI
#endif
  qmrhprint = .false.

  res = 0.1
  am = 0.05
  imass = 3
  iflag = 0

  anum(0) = 0.5
  do i = 1, ndiag
    anum(i) = real(exp(iunit*i*tau/ndiag))
    aden(i) = real(exp(-iunit*0.5*i*tau/ndiag))
  enddo

  ! dSdpi_ref is passed here simply to reduce the complexity of generate_starting_state
  call generate_starting_state(Phi, reqs_Phi, u, R, reqs_R, X2, reqs_X2, dSdpi_ref, Phi0_orig)

  ! call function
  do i = 1, timing_loops
    Phi0 = Phi0_orig
    call qmrherm(Phi, X2, res, itercg, am, imass, anum, aden, ndiag, iflag)
  end do

  ! check output
  if (generate) then
#ifdef MPI
    write_file(X2(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :), 'test_qmrherm_0_x.dat', MPI_Double_Complex)
    write_file(Phi0(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :), 'test_qmrherm_0_Phi0.dat', MPI_Double_Complex)
#else
    write (6, *) "Generation not possible"
    call exit(1)
#endif
  else
    read_file(X2_ref, 'test_qmrherm_0_x.dat', MPI_Double_Complex)
    read_file(Phi0_ref, 'test_qmrherm_0_Phi0.dat', MPI_Double_Complex)

    delta_x = X2(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - X2_ref
    delta_Phi0 = Phi0(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :) - Phi0_ref

    check_equality(itercg, 2, 'itercg', 'test_qmrherm_0')
    check_sum(delta_x, 50, 'X2', sum_delta_x, MPI_Double_Complex, 'test_qmrherm_0')
    check_max(delta_x, 0.5, 'X2', max_delta_x, MPI_Double_Precision, 'test_qmrherm_0')
    check_sum(delta_Phi0, 6e-12, 'Phi0', sum_delta_Phi0, MPI_Double_Complex, 'test_qmrherm_0')
    check_max(delta_Phi0, 1e-14, 'Phi0', max_delta_Phi0, MPI_Double_Precision, 'test_qmrherm_0')
  end if

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program test_qmrherm_0
