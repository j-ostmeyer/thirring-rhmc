#include "test_utils.fh"
program test_force
  use dwf3d_lib
  use trial
  use vector
  use dirac
  use gammamatrices
  use gforce
  use remez
  use remezg
  use counters
  use comms
  use comms4
  use comms5
  use comms6
  use qmrherm_module, only: phi0, qmrhprint => printall
  use test_utils
  implicit none

  ! general parameters
  logical :: generate = .false.
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  real :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)

  ! initialise function parameters
  complex(dp) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, 1)
  complex(dp) :: Phi0_orig(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, 25)
  real :: diff(ksizex_l, ksizey_l, ksizet_l, 3)
  complex(dp) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real :: sum_diff, max_diff

  integer :: imass, iflag, isweep, iter
  real :: res1, res2, am, u_variation
  real(dp) :: h, hg, hp, s, s_old

  integer :: i, j, l, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0

#ifdef MPI
  integer, dimension(12) :: reqs_R, reqs_X, reqs_U, reqs_Phi, reqs_Phi0
  integer :: ierr
  call init_MPI
#endif
  qmrhprint = .false.

  h = 0
  hg = 0
  hp = 0
  s = 0
  res1 = 1.0d-6
  res2 = 1.0d-8
  am = 0.05
  imass = 3
  iflag = 0
  isweep = 1
  iter = 0

  open (unit=36, file='remez2', status='old')
  open (unit=37, file='remez4', status='old')
  open (unit=38, file='remez2g', status='old')
  open (unit=39, file='remez4g', status='old')

  read (36, *) anum2(0)
  read (37, *) anum4(0)
  read (38, *) anum2g(0)
  read (39, *) anum4g(0)
  do i = 1, ndiag
    read (36, *) anum2(i), aden2(i)
    read (37, *) anum4(i), aden4(i)
  enddo
  do i = 1, ndiagg
    read (38, *) anum2g(i), aden2g(i)
    read (39, *) anum4g(i), aden4g(i)
  enddo
  read (36, *) bnum2(0)
  read (37, *) bnum4(0)
  read (38, *) bnum2g(0)
  read (39, *) bnum4g(0)
  do i = 1, ndiag
    read (36, *) bnum2(i), bden2(i)
    read (37, *) bnum4(i), bden4(i)
  enddo
  do i = 1, ndiagg
    read (38, *) bnum2g(i), bden2g(i)
    read (39, *) bnum4g(i), bden4g(i)
  enddo
!
  close (36)
  close (37)
  close (38)
  close (39)

  do j = 1, 4
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          do ithird = 1, kthird
            idx = ithird + (ip_x*ksizex_l + ix - 1)*kthird &
              & + (ip_y*ksizey_l + iy - 1)*kthird*ksize &
              & + (ip_t*ksizet_l + it - 1)*kthird*ksize*ksize &
              & + (j - 1)*kthird*ksize*ksize*ksizet
            Phi(ithird, ix, iy, it, j, 1) = 1.1*exp(iunit*idx*tau/idxmax)
            R(ithird, ix, iy, it, j) = 1.3*exp(iunit*idx*tau/idxmax)
            X(ithird, ix, iy, it, j) = 0.5*exp(1.0)*exp(iunit*idx*tau/idxmax)
            do l = 1, 25
              Phi0_orig(ithird, ix, iy, it, j, l) = 1.7*exp(1.0)*exp(iunit*idx*tau/idxmax) + l
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, R, 0, reqs_R)
  call start_halo_update_5(4, X, 0, reqs_X)
  call start_halo_update_6(4, 1, Phi, 1, reqs_Phi)
  call start_halo_update_6(4, 25, Phi0_orig, 2, reqs_Phi0)
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
          theta(ix, iy, it, j) = 1.9*exp(iunit*idx*tau/idxmax)
          pp(ix, iy, it, j) = -1.1*exp(iunit*idx*tau/idxmax)
          dSdpi(ix, iy, it, j) = tau*exp(iunit*idx*tau/idxmax)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(3, u, 3, reqs_u)
  call complete_halo_update(reqs_R)
  call complete_halo_update(reqs_X)
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_Phi0)
  call complete_halo_update(reqs_u)
#else
  call update_halo_6(4, 25, Phi0_orig)
  call update_halo_6(4, 1, Phi)
  call update_halo_5(4, R)
  call update_halo_5(4, X)
  call update_halo_4(3, u)
#endif

  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1

  call init_gammas()
  ! call function
  Phi0 = Phi0_orig
  h = 0
  hg = 0
  hp = 0
  s = 0
  max_qmr_iters = 10000
  !call gdbwait()
  call hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)
  s_old = s
  u_variation = 0.0001
  if (ip_global .eq. 1) then ! changing only one link in the whole lattice
    u(1, 1, 1, 1) = u(1, 1, 1, 1) + cmplx(0.0, u_variation/2)
    theta(1, 1, 1, 1) = theta(1, 1, 1, 1) + u_variation/2
  endif
  ! calculating force half way
  call force(Phi, res1, am, imass, isweep, iter)
  if (ip_global .eq. 1) then! changing only one link in the whole lattice
    u(1, 1, 1, 1) = u(1, 1, 1, 1) + cmplx(0.0, u_variation/2)
    theta(1, 1, 1, 1) = theta(1, 1, 1, 1) + u_variation/2
  endif
  call hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)

  check: block
    real(dp) :: action_difference
    real(dp) :: force
    if (ip_global .eq. 1) then
      action_difference = s - s_old
      force = dSdpi(1, 1, 1, 1)
      check_float_equality(force*u_variation, action_difference, 0.01, 'f*du', 'test_force')
    endif
  end block check

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
