#include "test_utils.fh"
program test_force
  use dwf3d_lib
  use trial
  use vector
  use dirac
  use gforce
  use remezg
  use avgitercounts
  use comms
  use qmrherm_module, only : phi0, qmrhprint => printall
  use test_utils
  implicit none

  ! general parameters
  logical :: generate = .false.
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8 * atan(1.0_8)

  real :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)

  ! initialise function parameters
  complex(dp) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, 1)
  complex(dp) :: Phi0_orig(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, 25)
  real :: diff(ksizex_l, ksizey_l, ksizet_l, 3)
  complex(dp) :: R(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  real :: sum_diff, max_diff

  integer :: imass, iflag, isweep, iter
  real :: res1,res2, am, u_variation
  real(dp) :: h, hg, hp, s, s_old

  integer :: i, j, l, ix, iy, it, ithird
  integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
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

  anum2g(0) = 0.5
  anum4g(0) = 0.51
  bnum2g(0) = 0.49
  bnum4g(0) = 0.53
  do i = 1, ndiagg
    anum2g(i) = 0.4 * real(exp(iunit * i * tau / ndiagg))
    aden2g(i) = 0.4 * real(exp(-iunit * 0.5 * i * tau / ndiagg))
    anum4g(i) = 0.41 * real(exp(iunit * i * tau / ndiagg))
    aden4g(i) = 0.41 * real(exp(-iunit * 0.5 * i * tau / ndiagg))
  enddo
  do j = 1,4
    do it = 1,ksizet_l
      do iy = 1,ksizey_l
        do ix = 1,ksizex_l
          do ithird = 1,kthird
            idx = ithird + (ip_x * ksizex_l + ix - 1) * kthird &
              & + (ip_y * ksizey_l + iy - 1) * kthird * ksize &
              & + (ip_t * ksizet_l + it - 1) * kthird * ksize * ksize &
              & + (j - 1) * kthird * ksize * ksize * ksizet
            Phi(ithird, ix, iy, it, j, 1) = 1.1 * exp(iunit * idx * tau / idxmax)
            R(ithird, ix, iy, it, j) = 1.3 * exp(iunit * idx * tau / idxmax)
            X(ithird, ix, iy, it, j) = 0.5 * exp(1.0) * exp(iunit * idx * tau / idxmax)
            do l = 1, 25
              Phi0_orig(ithird, ix, iy, it, j, l) = 1.7 * exp(1.0) * exp(iunit * idx * tau / idxmax) + l
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
  u = cmplx(1.0d0,2.0d-1)
  pp = -0.1
  dSdpi = 0
  theta = 2.0d-1
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

  call init(istart)
  ! call function
  Phi0 = Phi0_orig
  h = 0
  hg = 0
  hp = 0
  s = 0
  max_qmr_iters = 10000
  !call gdb_master_wait()
  call hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)
  print*,s,hg,hp,h
  s_old = s
  u_variation = 0.01
!  if(ip_global.eq.1)then
!    u(1,1,1,1) = u(1,1,1,1) + cmplx(0.0,u_variation/2)
!    theta(1,1,1,1) = theta(1,1,1,1) + u_variation/2
!  endif
  call force(Phi, res1, am, imass, isweep, iter)
  print*,s,hg,hp,h
!  if(ip_global.eq.1)then
!    u(1,1,1,1) = u(1,1,1,1) + cmplx(0.0,u_variation/2)
!    theta(1,1,1,1) = theta(1,1,1,1) + u_variation/2
!  endif
  call hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)
  print*,s,hg,hp,h

  if(ip_global.eq.1)then
    print*,"Action difference:",s-s_old
    print*,"Action:",s
    print*,"Old action",s_old
    ! check output
    print*,"Force*u_variation:",dSdpi(1,1,1,1)*u_variation
  endif

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
