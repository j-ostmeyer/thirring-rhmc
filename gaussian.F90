module gaussian
  use random
  use params

  real, parameter :: tpi = 2.0*acos(-1.0)

contains
  !**********************************************************************
  ! calculate vector of gaussian random numbers with unit variance
  ! to refresh momenta
  !   Numerical Recipes pp.203
  !**********************************************************************
#ifdef MPI
  subroutine gaussp(ps, reqs)
    use comms4, only: start_halo_update_4_real
#else
    subroutine gaussp(ps)
#endif
      real, intent(out) :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
#ifdef MPI
      integer, intent(out) :: reqs(12)
#endif
      integer ix, iy, it
      real :: theta
      !     write(6,1)
      !1   format(' Hi from gaussp')
      do it = 1, ksizet_l
        do iy = 1, ksizey_l
          do ix = 1, ksizex_l
            ps(ix, iy, it, 2) = sqrt(-2.0*log(rano(yran, idum, ix, iy, it)))
          end do
        end do
      end do
      do it = 1, ksizet_l
        do iy = 1, ksizey_l
          do ix = 1, ksizex_l
            theta = tpi*rano(yran, idum, ix, iy, it)
            ps(ix, iy, it, 1) = ps(ix, iy, it, 2)*sin(theta)
            ps(ix, iy, it, 2) = ps(ix, iy, it, 2)*cos(theta)
          end do
        end do
      end do
#ifdef MPI
      call start_halo_update_4_real(2, ps, 13, reqs)
#else
      call update_halo_4_real(2, ps)
#endif
      return
    end subroutine gaussp
    !**********************************************************************
    ! calculate vector of gaussian random numbers with unit variance
    ! to generate pseudofermion fields R
    !   Numerical Recipes pp.203
    !**********************************************************************
#ifdef MPI
    subroutine gauss0(ps, reqs)
      use comms4, only: start_halo_update_4_real
#else
      subroutine gauss0(ps)
#endif
        real, intent(out) :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
#ifdef MPI
        integer, intent(out) :: reqs(12)
#endif
        integer :: ix, iy, it
        real :: theta
        !     write(6,1)
        !1   format(' Hi from gauss0')
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              ps(ix, iy, it, 2) = sqrt(-log(rano(yran, idum, ix, iy, it)))
            end do
          end do
        end do
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              theta = tpi*rano(yran, idum, ix, iy, it)
              ps(ix, iy, it, 1) = ps(ix, iy, it, 2)*sin(theta)
              ps(ix, iy, it, 2) = ps(ix, iy, it, 2)*cos(theta)
            end do
          end do
        end do
#ifdef MPI
        call start_halo_update_4_real(2, ps, 14, reqs)
#else
        call update_halo_4_real(2, ps)
#endif
        return
      end subroutine gauss0

      end module gaussian
