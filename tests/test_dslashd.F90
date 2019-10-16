program test_dslashd
  use params
  use dwf3d_lib
  use dirac
  use comms
  use comms4
  use comms5
  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)
  complex(dp) :: acc_sum = 0.
  real(dp) :: acc_max = 0.

  ! common blocks to function

  ! initialise function parameters
  complex(dp) u(0:ksize + 1, 0:ksize + 1, 0:ksizet + 1, 3)
  complex(dp) Phi(kthird, 0:ksize + 1, 0:ksize + 1, 0:ksizet + 1, 4)
  complex(dp) Phiref(kthird, ksize, ksize, ksizet, 4)
  complex(dp) R(kthird, 0:ksize + 1, 0:ksize + 1, 0:ksizet + 1, 4)
  complex(dp) diff(kthird, ksize, ksize, ksizet, 4)

  real, parameter :: am = 0.05
  integer, parameter :: imass = 3

  integer :: i, j, l, ix, iy, it, ithird
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
#endif
  do j = 1, 3
    do it = 1, ksizet
      do iy = 1, ksize
        do ix = 1, ksize
          idx = ip_x*ksizex_l + ix &
            & + (ip_y*ksizey_l + iy - 1)*ksize &
            & + (ip_t*ksizet_l + it - 1)*ksize*ksize &
            & + (j - 1)*ksize*ksize*ksizet
          u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
        enddo
      enddo
    enddo
  enddo
  !      call update_halo_5(4, Phi)
#ifdef MPI
  call start_halo_update_4(3, u, 1, reqs_u)
  call complete_halo_update(reqs_R)
  call complete_halo_update(reqs_u)
#else
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
    call dslashd(Phi, R, u, am, imass)
#ifdef MPI
    call start_halo_update_5(4, Phi, 2, reqs_Phi)
    call complete_halo_update(reqs_Phi)
#else
    call update_halo_5(4, Phi)
#endif
  end do
  ! check output
  !      do i = 1,10
  !         j = 1 + i * (kvol - 1) / 10
  !         l = 1 + i * (4 - 1) / 10
  !         print *,'Phi(', j, ',', l, ') = ', Phi(j, l)
  !      enddo
  if (np_global .eq. 1) then
    open (3, file='test_dslashd.dat', form="unformatted", access="sequential")
    if (generate) then
      write (3) Phi(:, 1:ksize, 1:ksize, 1:ksizet, :)
    else
      read (3) Phiref

      diff = Phi(:, 1:ksize, 1:ksize, 1:ksizet, :) - Phiref
      print *, 'sum delta = ', sum(diff)
      print *, 'max delta = ', maxval(abs(diff))
    end if
  end if
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
