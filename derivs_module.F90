module derivs_module
  use mpi
  use comms_common
  implicit none

contains

  subroutine derivs(R, X2, anum, iflag, am, imass)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, dp
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag
    integer, intent(in) :: imass
    real, intent(in) :: am

#ifdef GENERATE_WITH_SHAMIR
    call derivs_shamir(R, X2, anum, iflag, am, imass)
#endif

#ifdef GENERATE_WITH_WILSON
    call derivs_wilson(R, X2, anum, iflag, am, imass)
#endif

  end subroutine derivs 

  subroutine derivs_shamir(R, X2, anum, iflag, am, imass)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, dp
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag
    integer, intent(in) :: imass
    real, intent(in) :: am

    call derivs_shared(R, X2, anum, iflag, am, imass)
  end subroutine derivs_shamir

  subroutine derivs_wilson(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, dp
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag
    integer, intent(in) :: imass
    real, intent(in) :: am

    complex(dp) :: cmult
    integer :: il, ierr
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real sumdS_S,sumdS_W,absdS_S,absdS_W

    call derivs_shared(R, X2, anum, iflag, am, imass)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This Calculates a sum for the shamir version of derivs before any wilson aspects have been calculated
    ! Safe to comment out
    sumdS_S=sum(dSdpi)
    absdS_S=sum(abs(dSdpi))
    call MPI_AllReduce(MPI_In_Place, sumdS_S, 1, MPI_Real, MPI_Sum, comm, ierr)
    call MPI_AllReduce(MPI_In_Place, absdS_S, 1, MPI_Real, MPI_Sum, comm, ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call derivsWilsonMass(R, X2, anum, iflag, am, imass)

    cmult=cmplx(1.0,0.0)

    do il=1,kthird_l-1
      if (iflag.eq.1) then
        call derivsDMTest(R, X2, anum, iflag, am, imass,il,il+1,cmult)
        call derivsDPTest(R, X2, anum, iflag, am, imass,il+1,il,cmult)
      elseif (iflag.eq.0) then
        call derivsPDTest(R, X2, anum, iflag, am, imass,il,il+1,cmult)
        call derivsMDTest(R, X2, anum, iflag, am, imass,il+1,il,cmult)
      endif
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Safe to comment out
    sumdS_W=sum(dSdpi)
    absdS_W=sum(abs(dSdpi))
    call MPI_AllReduce(MPI_In_Place, sumdS_W, 1, MPI_Real, MPI_Sum, comm, ierr)
    call MPI_AllReduce(MPI_In_Place, absdS_W, 1, MPI_Real, MPI_Sum, comm, ierr)

    if (ip_global .eq. 0) then
        open(unit=105,file="fort.105",position="append")
        write(105,'(5f20.12)') sumdS_S,sumdS_W,absdS_S,absdS_W,absdS_W/absdS_S
        close(105)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine derivs_wilson

  subroutine derivs_shared(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    !      complex, intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    !      complex, intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp) :: dSdpi_tmp(ksizex_l, ksizey_l, ksizet_l, 3)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer :: ierr

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
                        * sum(aimag(conjg(R(1:kthird_l, ix, iy, it, idirac)) &
                                        * X2(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                            - aimag(conjg(R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                        * X2(1:kthird_l, ix, iy, it, idirac)))
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
                            * sum(aimag(gamval(mu, idirac) &
                                * (conjg(R(1:kthird_l, ix, iy, it, idirac)) &
                                  * X2(1:kthird_l, ix + ixup, iy + iyup, it + itup, igork1) &
                                 + conjg(R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                  * X2(1:kthird_l, ix, iy, it, igork1))))
              enddo
            enddo
          enddo

        else

          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                dSdpi_tmp(ix, iy, it, mu) = dSdpi_tmp(ix, iy, it, mu) &
                            - tzi_real &
                            * sum(aimag(gamval(mu, idirac) &
                                * (conjg(R(1:kthird_l, ix, iy, it, idirac)) &
                                  * X2(1:kthird_l, ix + ixup, iy + iyup, it + itup, igork1) &
                                 + conjg(R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac)) &
                                  * X2(1:kthird_l, ix, iy, it, igork1))))
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
  end subroutine derivs_shared

  subroutine derivsWilsonMass(R, X2, anum, iflag, am, imass)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, dp
    use comms
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass

    complex(dp) :: cmult
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer idxL,idxR


    ! upper right mass 
    if (imass.eq.1) then
      cmult=cmplx(-am,0)
    elseif (imass.eq.3) then
      cmult=cmplx(0,-am)
    endif

    if (iflag.eq.1) then
      idxL=1
      idxR=kthird_l
      call derivsDPTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    elseif (iflag.eq.0) then
      idxL=1
      idxR=kthird_l
      call derivsMDTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    end if

    ! lower left mass 
    if (imass.eq.1) then
      cmult=cmplx(-am,0)
    elseif (imass.eq.3) then
      cmult=cmplx(0,am)
    endif

    if (iflag.eq.1) then
      idxL=kthird_l
      idxR=1
      call derivsDMTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    elseif (iflag.eq.0) then
      idxL=kthird_l
      idxR=1
      call derivsPDTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    endif

    return
  end subroutine derivsWilsonMass


  subroutine derivsDPTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use comms
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR
    complex(dp),intent(in) :: cmult

    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    sliceL(:,:,:,:)=R(idxL,:,:,:,:)
    sliceR=cmplx(0.0,0.0)
    sliceR(:,:,:,1:2)=X2(idxR,:,:,:,1:2)
    call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

    return
  end subroutine derivsDPTest


  subroutine derivsDMTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use comms
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR
    complex(dp),intent(in) :: cmult

    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    sliceL(:,:,:,:)=R(idxL,:,:,:,:)
    sliceR=cmplx(0.0,0.0)
    sliceR(:,:,:,3:4)=X2(idxR,:,:,:,3:4)
    call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

    return
  end subroutine derivsDMTest


  subroutine derivsPDTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use comms
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR
    complex(dp),intent(in) :: cmult

    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    sliceL=cmplx(0.0,0.0)
    sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)
    sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
    call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

    return
  end subroutine derivsPDTest


  subroutine derivsMDTest(R, X2, anum, iflag, am, imass,idxL,idxR,cmult)
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use comms
    implicit none
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR
    complex(dp),intent(in) :: cmult

    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    sliceL=cmplx(0.0,0.0)
    sliceL(:,:,:,3:4)=R(idxL,:,:,:,3:4)
    sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
    call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

    return
  end subroutine derivsMDTest


  subroutine derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird_l, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use comms
    implicit none
    complex(dp), intent(in) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass
    complex(dp), intent(in) :: cmult

    real :: tzi_real
    integer :: ixup,iyup,itup,ix, iy, it, idirac, mu
    integer :: igork1
    complex(dp) :: expp, expm
    real :: dmult

    tzi_real = 2*real(anum)
    expp=1d0
    expm=1d0

    do mu = 1, 3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4 
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
                dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + tzi_real*real(akappa)*( &
              &   aimag( cmult* conjg(sliceL(ix, iy, it, idirac))* expp * &
              &           sliceR(ix + ixup, iy + iyup, it + itup, idirac) ) &
              & - aimag( cmult* conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
              &           sliceR(ix, iy, it, idirac) ) )
            enddo
          enddo
        enddo

        dmult=1.0
        if (iflag.eq.1) then
          dmult=-1.0
        endif

        igork1 = gamin(mu, idirac)
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
                dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + tzi_real*dmult*aimag( cmult*gamval(mu, idirac)* &
              &  ( conjg(sliceL(ix, iy, it, idirac))* expp * &
              &        sliceR(ix + ixup, iy + iyup, it + itup, igork1) &
              &  + conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
              &        sliceR(ix, iy, it, igork1)) )
            enddo
          enddo
        enddo
      enddo ! idirac
    enddo ! mu

    return
  end subroutine derivSlice


end module derivs_module
