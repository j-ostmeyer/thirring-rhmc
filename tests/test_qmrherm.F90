#include "test_utils.fh"
program test_qmrherm_0
  use qmrherm_module, only: qmrherm, phi0, qmrhprint => printall
  use trial, only: u
  use vector, only: X
  use gammamatrices
  use gforce
  use params
  use comms
  use comms4
  use comms5
  use comms6
  use test_utils
  implicit none

  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)
  integer :: imass, iflag, i, itercg, ierr
  character(len=4) :: iflag_char
  character(len=14) :: test_name
  real :: res, am
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: Phi0_ref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, 25)
  complex(dp) :: delta_Phi0(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, 25)
  complex(dp) :: Phi0ToWrite(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, 25)
  complex(dp) :: Phi0_orig(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, 25)
  complex(dp) :: xToWrite(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: x_ref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: delta_x(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: sum_delta_x, sum_delta_Phi0
  real(dp) :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)
  real(dp) :: anum(0:ndiag), aden(ndiag)
  real(dp) :: max_delta_x, max_delta_Phi0
  integer, dimension(16) :: reqs_R, reqs_Phi, reqs_Phi0, reqs_X
  integer, dimension(12) :: reqs_u

#ifdef MPI
  call init_MPI
#endif

  qmrhprint = .false.

  res = 0.1
  am = 0.05
  imass = 3

  do iflag = 0, 3
    ! Setup iflag_char
    iflag_char = ''
    write(iflag_char, '(I1)') iflag
    if (ip_global == 0) then
      print *, ' iflag: ', iflag_char
    end if

    test_name = "test_qmrherm_"//iflag_char

    anum(0) = 0.5
    do i = 1, ndiag
      anum(i) = real(exp(iunit*i*tau/ndiag))
      aden(i) = real(exp(-iunit*0.5*i*tau/ndiag))
    enddo

    call generate_starting_state(Phi, reqs_Phi, u, R, reqs_R, X, reqs_X, dSdpi_ref, Phi0_orig)

    ! call function
    Phi0 = Phi0_orig
    call qmrherm(Phi, X, res, itercg, am, imass, anum, aden, ndiag, iflag)

    ! check output
    if (generate) then
#ifdef MPI
      xToWrite(:, :, :, :, :) = x(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
      Phi0ToWrite(:, :, :, :, :, :) = Phi0(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :)
      write_file(xToWrite(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :), test_name//'_x.dat', MPI_Double_Complex)
      write_file(Phi0ToWrite(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :), test_name//'_Phi0.dat', MPI_Double_Complex)
#else
      write (6, *) "Generation not possible"
      call exit(1)
#endif
    else
      read_file(x_ref, test_name//'_x.dat', MPI_Double_Complex)
      read_file(Phi0_ref, test_name//'_Phi0.dat', MPI_Double_Complex)

      delta_x = x(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - x_ref
      delta_Phi0 = Phi0(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, :) - Phi0_ref

      check_equality(itercg, 2, 'itercg', test_name)
      check_sum(delta_x, 50, 'x', sum_delta_x, MPI_Double_Complex, test_name)
      check_max(delta_x, 0.5, 'x', max_delta_x, MPI_Double_Precision, test_name)
      check_sum(delta_Phi0, 6e-12, 'Phi0', sum_delta_Phi0, MPI_Double_Complex, test_name)
      check_max(delta_Phi0, 1e-14, 'Phi0', max_delta_Phi0, MPI_Double_Precision, test_name)
    end if
  end do

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program test_qmrherm_0
