#include "test_utils.fh"
program test_qmrherm_4
      use inverter_checks
      use qmrherm_scratch
      use dwf3d_lib
      use trial
      use phizero
      use gforce
      use comms
      use test_utils
      implicit none

! general parameters
      complex, parameter :: iunit = cmplx(0, 1)
      real(dp), parameter :: tau = 8 * atan(1.0_8)

! common blocks to function
      

! initialise function parameters
      complex(dp) Phi(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), allocatable :: Phi0_ref(:, :, :, :, :, :)
      complex(dp), allocatable :: Phi0_orig(:, :, :, :, :, :)
      complex(dp), allocatable :: delta_Phi0(:, :, :, :, :, :)
      complex(dp) :: xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      !complex(dp) :: R(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      integer :: imass, iflag, isweep, iter
      real(dp) :: anum(0:ndiag), aden(ndiag), diffnorm2
      real :: res, am, anumf, adenf
      integer :: itercg
      
      integer :: i, j, l, ix, iy, it, ithird
      integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
      integer :: idx = 0, idiag

#ifdef MPI
      type(MPI_Request), dimension(12) :: reqs_R, reqs_U, reqs_Phi, reqs_Phi0
      type(MPI_Request), dimension(12) :: reqs_xin
      call init_MPI
#endif

      allocate(Phi0_ref(kthird, ksizex_l, ksizey_l, ksizet_l, 4, 25))
      allocate(Phi0_orig(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, 25))
      allocate(delta_Phi0(kthird, ksizex_l, ksizey_l, ksizet_l, 4, 25))

      res = 1e-4
      am = 0.05
      imass = 3
      iflag = 0
      isweep = 1
      iter = 0

      open(unit=36,file='../remez2',status='old')
      read(36,*) anum(0)
      do i=1,ndiag
          read(36,*) anum(i),aden(i)
      enddo
      close(36)
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
                     R(ithird, ix, iy, it, j) = 1.3 * exp(iunit * idx * tau / idxmax)
                     do l = 1, 25
                        Phi0_orig(ithird, ix, iy, it, j, l) = &
                             & 1.7 * exp(1.0) * exp(iunit * idx * tau / idxmax) + l
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
      do j = 1,3
         do it = 1,ksizet_l
            do iy = 1,ksizey_l
               do ix = 1,ksizex_l
                  idx = ip_x * ksizex_l + ix - 1 &
                       & + (ip_y * ksizey_l + iy - 1) * ksize &
                       & + (ip_t * ksizet_l + it - 1) * ksize * ksize &
                       & + (j - 1) * ksize * ksize * ksizet
                  u(ix, iy, it, j) = exp(iunit * idx * tau / idxmax)
                  dSdpi(ix, iy, it, j) = real(tau * exp(iunit * idx * tau / idxmax), sp)
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
      call update_halo_6(4, 25, Phi0_orig)
      call update_halo_5(4, Phi)
      call update_halo_5(4, R)
      call update_halo_4(3, u)
#endif
! initialise common variables
      beta = 0.4
      am3 = 1.0
      ibound = -1
      
      call init(istart)
! call function
      Phi0 = Phi0_orig
      call qmrherm(Phi, res, itercg, am, imass, anum, aden, ndiag, iflag, isweep, iter)
! check output
      do idiag=1,ndiag
          xin(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(:,:,:,:,:,idiag)
#ifdef MPI
          call start_halo_update_5(4, xin, 10, reqs_xin)
          call complete_halo_update(reqs_xin)
#else
          call update_halo_5(4, xin)
#endif
          adenf = aden(idiag)
          anumf = anum(idiag)
          call dirac_op_shifted(xout,xin,am,imass,adenf,anumf)
          call check_diff(diffnorm2,xout,Phi) 
          write(6,*) "Diffnorm2 - shift " ,idiag,"=", diffnorm2


      enddo
#ifdef MPI
      call MPI_Finalize
#endif
end program test_qmrherm_4
