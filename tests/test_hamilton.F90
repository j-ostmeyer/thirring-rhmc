program test_hamilton
  use dwf3d_lib
  use counters
  use trial
  use qmrherm_module, only: phi0, qmrhprint => printall
  use dirac
  use gammamatrices
  use gforce
  use remez
  use remezg
  use params
  use comms
  use comms4
  use comms5
  use comms6
  implicit none

  ! general parameters
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  complex(dp) Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp), allocatable :: Phi0_ref(:, :, :, :, :, :)
  complex(dp), allocatable :: Phi0_orig(:, :, :, :, :, :)
  complex(dp) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

  integer :: imass, iflag, isweep, iter
  real :: res2, am
  real(dp) :: h, hg, hp, s

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

  hg = 0.
  hp = 0.
  h = 0.
  s = 0.

  res2 = 0.1
  am = 0.05
  imass = 3
  iflag = 0
  isweep = 1
  iter = 0

  anum2(0) = 0.5
  anum4(0) = 0.51
  bnum2(0) = 0.49
  bnum4(0) = 0.53
  do i = 1, ndiag
    anum2(i) = 0.4*exp(iunit*i*tau/ndiag)
    aden2(i) = 0.4*exp(-iunit*0.5*i*tau/ndiag)
    anum4(i) = 0.41*exp(iunit*i*tau/ndiag)
    aden4(i) = 0.41*exp(-iunit*0.5*i*tau/ndiag)
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
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_Phi0)
  call complete_halo_update(reqs_u)
#else
  call update_halo_6(4, 25, Phi0)
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
    h = 0
    hg = 0
    hp = 0
    s = 0
    call hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)
  end do
  ! check output
  if (ip_global .eq. 0) then
    if (ancghpv .ne. 2.0) then
      print *, "ancghpv looks wrong:", ancghpv
    end if
    if (ancgh .ne. 2.0) then
      print *, "ancgh looks wrong:", ancgh
    end if
    if (hg .lt. 3695.10 .or. hg .gt. 3695.11) then
      print *, "hg looks wrong:", hg
    end if
    if (hp .lt. 3096.31 .or. hp .gt. 3096.32) then
      print *, "hp looks wrong:", hp
    end if
    ! Changes to summation in calculation of alphatild inside qmrherm mean that this
    ! test gives a slightly different value to in the original code, hence the looser check.
    ! In principle using sum() rather than direct iterative summation gives a more precise answer
    if (s .lt. 144000. .or. s .gt. 145000.) then
      print *, "s looks wrong:", s
    end if
  end if
  call MPI_Finalize(ierr)
end program
