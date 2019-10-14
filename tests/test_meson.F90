#include "test_utils.fh"
program test_measure
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

  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real(dp), parameter :: tau = 8*atan(1.0_8)
  complex(dp) :: cferm1(0:ksizet - 1), cferm2(0:ksizet - 1)
  real(dp) :: cpm(0:ksizet - 1), cmm(0:ksizet - 1)
  complex(dp) :: cferm1_ref(0:ksizet - 1), cferm2_ref(0:ksizet - 1)
  real(dp) :: cpm_ref(0:ksizet - 1), cmm_ref(0:ksizet - 1), maxreldiff

  ! initialise function parameters
  real :: aviter
  integer :: iflag = 0

  complex(dp) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  integer :: i, j, ix, iy, it, ithird, idx
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  real :: res, am
  integer :: imass, isweep, itercg, iter
#ifdef MPI
  integer, dimension(12) :: reqs_Phi, reqs_u
  integer :: ierr

  call init_MPI
#endif
  seed = 4139764973254.0
  idum = -1
  call rranset(seed, 1, 1, 1)
  seed = rano(yran, idum, 1, 1, 1)

  res = 0.1
  am = 0.05
  imass = 3
  iflag = 0
  isweep = 1
  iter = 0

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
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, Phi, 0, reqs_Phi)
#endif
  idx = 0
  do j = 1, 3
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          idx = ip_x*ksizex_l + ix &
            & + (ip_y*ksizey_l + iy - 1)*ksize &
            & + (ip_t*ksizet_l + it - 1)*ksize*ksize &
            & + (j - 1)*ksize*ksize*ksizet
          u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(3, u, 1, reqs_u)
  call complete_halo_update(reqs_Phi)
  call complete_halo_update(reqs_u)
#else
  call update_halo_5(4, Phi)
  call update_halo_4(3, u)
#endif

  ! initialise common variables
  beta = 0.4
  am3 = 1.0
  ibound = -1
  seed = rano(yran, idum, 1, 1, 1)
  call init_gammas()
  seed = rano(yran, idum, 1, 1, 1)

  ! call function
  do i = 1, timing_loops
    x = (0.D0, 0.D0)
    call meson(cpm, cmm, cferm1, cferm2, res, itercg, aviter, am, imass)
  end do

  if (generate) then
    open (3, file='test_meson.dat', form="unformatted", access="sequential")
#ifdef MPI
    if (ip_global .eq. 0) then
#endif
      write (3) cferm1, cferm2, cpm, cmm
      print *, "cferm1"
      print *, cferm1
      print *, "cferm2"
      print *, cferm2
      print *, "cpm"
      print *, cpm
      print *, "cmm"
      print *, cmm
#ifdef MPI
    endif! if(ip_global .eq. 0) then
#endif
  else
    open (3, file='test_meson.dat', form="unformatted", access="sequential", status='old')
#ifdef MPI
    if (ip_global .eq. 0) then
#endif
      read (3) cferm1_ref, cferm2_ref, cpm_ref, cmm_ref
      maxreldiff = maxval(2*abs(cferm1_ref - cferm1)/abs(cferm1_ref + cferm1))
      if (maxreldiff .gt. 1.0e-07) then
        print *, 'maxval(2*abs(cferm1_ref - cferm1)/abs(cferm1_ref + cferm1))=', maxreldiff
      endif
      maxreldiff = maxval(2*abs(cferm2_ref - cferm2)/abs(cferm2_ref + cferm2))
      if (maxreldiff .gt. 1.0e-06) then
        print *, 'maxval(2*abs(cferm2_ref - cferm2)/abs(cferm2_ref + cferm2))=', maxreldiff
      endif
      maxreldiff = maxval(2*abs(cpm_ref - cpm)/abs(cpm_ref + cpm))
      if (maxreldiff .gt. 1.0e-07) then
        print *, 'maxval(2*abs(cpm_ref - cpm)/abs(cpm_ref + cpm))=', maxreldiff
      endif
      maxreldiff = maxval(2*abs(cmm_ref - cmm)/abs(cmm_ref + cmm))
      if (maxreldiff .gt. 1.0e-07) then
        print *, 'maxval(2*abs(cmm_ref - cmm)/abs(cmm_ref + cmm))=', maxreldiff
      endif

#ifdef MPI
    endif! if(ip_global .eq. 0) then
    call MPI_Finalize(ierr)
#endif
  end if
  close (3)
end program
