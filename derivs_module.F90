module derivs_module
  use mpi
  use comms_common, only: comm_grp_third
  implicit none

contains

  subroutine derivs(R, X2, anum, iflag)
    use gforce, only: dSdpi
    use dirac, only: kdelta, gamval, gamin
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    !      complex, intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    !      complex, intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag
    integer :: ierr
    real(dp)  :: dSdpi_tmp(ksizex_l, ksizey_l, ksizet_l, 3)

    !      complex(dp) :: tzi
    real :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, mu
    integer :: igork1

    !     write(6,111)
    !111 format(' Hi from derivs')

    !     dSdpi_tmp=dSdpi_tmp-Re(Rdagger *(d(Mdagger)dp)* X2)
    !     Cf. Montvay & Muenster (7.215)
    !      tzi=cmplx(0.0,2*anum)
    tzi_real = 2*real(anum)
    !     factor of 2 picks up second term in M&M (7.215)
    dSdpi_tmp = 0.0

    do mu = 1, 3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac = 1, 4
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              dSdpi_tmp(ix, iy, it, mu) = dSdpi_tmp(ix, iy, it, mu) &
                                          + tzi_real*real(akappa) &
                                          *sum(aimag(conjg(R(1:kthird_l, ix, iy, it, idirac)) &
                                                     *X2(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                               - aimag(conjg(R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                                       *X2(1:kthird_l, ix, iy, it, idirac)))
            enddo
          enddo
        enddo

        igork1 = gamin(mu, idirac)

        if (iflag .eq. 0) then

          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                dSdpi_tmp(ix, iy, it, mu) = dSdpi_tmp(ix, iy, it, mu) &
                                            + tzi_real &
                                            *sum(aimag(gamval(mu, idirac) &
                                                       *(conjg(R(1:kthird_l, ix, iy, it, idirac)) &
                                                         *X2(1:kthird_l, ix + ixup, iy + iyup, it + itup, igork1) &
                                                         + conjg(R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                                         *X2(1:kthird_l, ix, iy, it, igork1))))
              enddo
            enddo
          enddo

        else

          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                dSdpi_tmp(ix, iy, it, mu) = dSdpi_tmp(ix, iy, it, mu) &
                                            - tzi_real &
                                            *sum(aimag(gamval(mu, idirac) &
                                                       *(conjg(R(1:kthird_l, ix, iy, it, idirac)) &
                                                         *X2(1:kthird_l, ix + ixup, iy + iyup, it + itup, igork1) &
                                                         + conjg(R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                                         *X2(1:kthird_l, ix, iy, it, igork1))))
              enddo
            enddo
          enddo

        endif

      enddo
    enddo

#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, dSdpi_tmp, ksizex_l*ksizey_l*ksizet_l*3, &
                       MPI_DOUBLE_PRECISION, MPI_Sum, comm_grp_third, ierr)
#endif
    dSdpi = dSdpi + dSdpi_tmp

    return
  end subroutine derivs

end module derivs_module
