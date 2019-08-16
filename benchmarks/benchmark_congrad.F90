program benchmark_congrad
  use dwf3d_lib
  use trial, only: u
  use vector, only : x
  use dirac
  use comms
  use comms4
  use comms5
  use comms6
  use measure_module
  implicit none

  ! general parameters
  complex, parameter :: iunit = cmplx(0, 1)
  real*8, parameter :: tau = 8 * atan(1.0_8)



  ! initialise function parameters
  complex(dp) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: Phi0(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: X0(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

  integer :: imass, iflag, isweep, iter
  real :: res, am
  integer :: itercg

  integer :: i, j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
  integer :: idx = 0
  integer :: iterations
  !integer :: t1i(2), t2i(2)
  double precision :: t1i, t2i
  double precision :: dt
#ifdef MPI
  integer, dimension(12) :: reqs_X, reqs_Phi, reqs_u
  integer :: ierr

  call init_MPI
#endif

  res = 1e-14
  am = 0.05
  imass = 3
  iflag = 0
  isweep = 1
  iter = 0

  do j = 1,4
    do it = 1,ksizet_l
      do iy = 1,ksizey_l
        do ix = 1,ksizex_l
          do ithird = 1,kthird
            idx = ithird + (ip_x * ksizex_l + ix - 1) * kthird &
              & + (ip_y * ksizey_l + iy - 1) * kthird * ksize &
              & + (ip_t * ksizet_l + it - 1) * kthird * ksize * ksize &
              & + (j - 1) * kthird * ksize * ksize * ksizet
            Phi(ithird, ix, iy, it, j) = 1.1 * exp(iunit * idx * tau / idxmax)
            X(ithird, ix, iy, it, j) = 0.5 * exp(1.0) * exp(iunit*idx*tau/idxmax)
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, X, 0, reqs_X)
  call start_halo_update_5(4, Phi, 0, reqs_Phi)
#endif
  idx = 0
  do j = 1,3
    do it = 1,ksizet_l
      do iy = 1,ksizey_l
        do ix = 1,ksizex_l
          idx = ip_x * ksizex_l + ix &
            & + (ip_y * ksizey_l + iy - 1) * ksize &
            & + (ip_t * ksizet_l + it - 1) * ksize * ksize &
            & + (j - 1) * ksize * ksize * ksizet
          u(ix, iy, it, j) = exp(iunit * idx * tau / idxmax)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(3, u, 1, reqs_u)
  call complete_halo_update(reqs_X)
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_u)
#else
  call update_halo_5(4, Phi)
  call update_halo_5(4, X)
  call update_halo_4(3, u)
#endif
  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1

  call init(istart)
  ! call function
  Phi0 = Phi
  X0 = X

  !call gettimeofday(t1i,ierr)
  t1i = MPI_Wtime()
  do i = 1,timing_loops
    Phi=Phi0
    X = X0
    call congrad(Phi, res, itercg, am, imass,iterations)
#ifdef MPI
    if(ip_global.eq.0) then
#endif
      print *, "Congrad iterations, res:", iterations,res
#ifdef MPI
    endif
#endif
  end do
  !call gettimeofday(t2i,ierr)
  t2i = MPI_Wtime()
  !dt = t2i(1)-t1i(1) + 1.0d-6 * (t2i(2)-t1i(2))
  dt = t2i-t1i
  dt = dt / timing_loops / iterations
#ifdef MPI
  if(ip_global.eq.0) then
#endif
    print *, "Time per iteration:" , dt
#ifdef MPI
  endif
  call MPI_Finalize(ierr)
#endif

end program benchmark_congrad
