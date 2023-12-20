#include "test_utils.fh"
program test_measure
  use gdbhook
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
  use gdbhook
  implicit none

  integer :: imass_index, imass, ierr, itercg
  character(len=4) :: imass_char
  real :: res, am, aviter
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  integer, dimension(16) :: reqs_Phi
  
  ! expected test values
  real, dimension(3) :: expected_aviter = (/3,6,6/)

#ifdef MPI
  call init_MPI
  call gdb_wait()
#endif

  seed = 4139764973254.0
  idum = -1
  call rranset(seed, 1, 1, 1)
  seed = rano(yran, idum, 1, 1, 1)

  res = 0.1
  am = 0.05
  imass = 3

  do imass_index = 1, size(imasses)
    imass_char = ''
    imass = imasses(imass_index)
    write(imass_char, '(I1)') imass
    if (ip_global == 0) then
      print *, ' imass: ', imass_char
    end if
#ifdef MPI
    call MPI_Barrier(comm, ierr)
#endif
    call generate_starting_state(Phi, reqs_Phi, u)

    ! call function
    x = (0.D0, 0.D0)
    call meson(res, itercg, aviter, am, imass)

    check_equality(itercg, 3, 'itercg', 'test_meson')
    check_equality(aviter, expected_aviter(imass_index), 'aviter', 'test_meson')
  end do

#ifdef MPI
    call MPI_Finalize(ierr)
#endif
  close (3)
end program
