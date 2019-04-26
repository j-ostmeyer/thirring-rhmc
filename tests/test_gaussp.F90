#include "test_utils.fh"
program test_gaussp
  use dwf3d_lib
  use gauge
  use gaussian
  use comms
  use random
  implicit none

  real :: ps(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
  integer :: ix, iy, it, ix2, iy2, it2, i, j
  integer, dimension(4) :: duplicate_position1, duplicate_position2
  logical :: has_duplicates = .false.
  real :: sumps, maxps, minps

  ! initialise MPI
#ifdef MPI
  integer :: reqs(12)
  integer :: ierr
  call init_MPI
#endif

  seed = 1.0
  call init_random(seed)

  ! call function
#ifdef MPI
  call gaussp(ps, reqs)
  call complete_halo_update(reqs)
#else
  call gaussp(ps)
#endif

  ! check output
  sumps = sum(ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
  maxps = maxval(ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
  minps = minval(ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
  outer: do it = 1, ksizet_l
    do it2 = 1, ksizet_l
      do iy = 1, ksizey_l
        do iy2 = 1, ksizey_l
          do ix = 1, ksizex_l
            do ix2 = 1, ksizex_l
              do i = 1, 2
                do j = 1, 2
                  if (ix .eq. ix2 .and. iy .eq. iy2 .and. it .eq. it2 &
                    .and. i .eq. j) cycle
                  if (ps(ix, iy, it, i) .eq. ps(ix2, iy2, it2, j)) then
                    has_duplicates = .true.
                    duplicate_position1 = (/ ix, iy, it, i /)
                    duplicate_position2 = (/ ix2, iy2, it2, j /)
                    exit outer
                  endif
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do outer

#ifdef MPI
  call MPI_AllReduce(MPI_IN_PLACE, sumps, 1, MPI_REAL, MPI_SUM, comm,ierr)
  call MPI_AllReduce(MPI_IN_PLACE, maxps, 1, MPI_REAL, MPI_MAX, comm,ierr)
  call MPI_AllReduce(MPI_IN_PLACE, minps, 1, MPI_REAL, MPI_MIN, comm,ierr)
  call MPI_AllReduce(MPI_IN_PLACE, has_duplicates, 1, MPI_LOGICAL, MPI_LOR, comm,ierr)
#endif
#ifndef SITE_RANDOM
#if defined(MPI) && NP_T > 1 && NP_X > 1 && NP_Y > 1
  if (ip_global .eq. 0) then
    print *, 'Unable to check results in parallel without using site_random'
  end if
#else
  check_float_equality(sumps, -22.13475, 0.001, 'sum', 'test_gaussp')
  check_float_equality(maxps, 3.914956, 0.001, 'max', 'test_gaussp')
  check_float_equality(minps, -3.558903, 0.001, 'min', 'test_gaussp')
#endif
#else
  check_float_equality(sumps, 20.4506, 0.001, 'sum', 'test_gaussp')
  check_float_equality(maxps, 3.67040, 0.001, 'max', 'test_gaussp')
  check_float_equality(minps, -3.452905, 0.001, 'min', 'test_gaussp')
  if (ip_global .eq. 0) then
    if (has_duplicates) then
      print *, 'duplicate random numbers observed at:'
      write(*, '(4i3)') duplicate_position1
      write(*, '(4i3)') duplicate_position2
    end if
  end if
#endif
  if (ip_global .eq. 0) then
    if (has_duplicates) then
      print *, 'duplicate random numbers observed at:'
      write(*, '(4i3)') duplicate_position1
      write(*, '(4i3)') duplicate_position2
    end if
  end if
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
