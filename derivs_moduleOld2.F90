#include "kernel.h"

module derivs_module
  implicit none

contains

  subroutine derivs(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
!    logical,parameter :: COMPACT=.true.
    !      complex, intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    !      complex, intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass

    !      complex(dp) :: tzi
    real :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, mu
    integer :: igork1
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer idxL,idxR
    real sumdS_S,sumdS_W,absdS_S,absdS_W
    integer ierr
    !
    !     write(6,111)
    !111 format(' Hi from derivs')

    !     dSdpi=dSdpi-Re(Rdagger *(d(Mdagger)dp)* X2)
    !     Cf. Montvay & Muenster (7.215)
    !      tzi=cmplx(0.0,2*anum)
    tzi_real = 2*real(anum)
    !     factor of 2 picks up second term in M&M (7.215)
    !
    expp=1d0
    expm=1d0
    mult=-1d0
    do mu = 1, 3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac = 1, 4
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              if (COMPACT) then
                expp=exp(cmplx(0,theta(ix,iy,it,mu)))
                expm=exp(mult*cmplx(0,theta(ix,iy,it,mu)))
              endif
              dSdpi(ix, iy, it, mu) = &
               &     dSdpi(ix, iy, it, mu) + tzi_real*real(akappa)*sum(aimag( &
                &       conjg(R(:, ix, iy, it, idirac))* expp * &
                &         X2(:, ix + ixup, iy + iyup, it + itup, idirac)) &
          &     - aimag(conjg(R(:, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &         X2(:, ix, iy, it, idirac)))
            enddo
          enddo
        enddo
        !
        igork1 = gamin(mu, idirac)
        if (iflag .eq. 0) then ! this is dagger!!!
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                if (COMPACT) then
                  expp=exp(cmplx(0,theta(ix,iy,it,mu)))
                  expm=exp(mult*cmplx(0,theta(ix,iy,it,mu)))
                endif
                dSdpi(ix, iy, it, mu) = &
         &      dSdpi(ix, iy, it, mu) + tzi_real*sum(aimag(gamval(mu, idirac)* &
                  &(conjg(R(:, ix, iy, it, idirac))* expp * &
                  &        X2(:, ix + ixup, iy + iyup, it + itup, igork1) &
                  &+ conjg(R(:, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                  &             X2(:, ix, iy, it, igork1))))
              enddo
            enddo
          enddo
        else
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                if (COMPACT) then
                  expp=exp(cmplx(0,theta(ix,iy,it,mu)))
                  expm=exp(mult*cmplx(0,theta(ix,iy,it,mu)))
                endif
                dSdpi(ix, iy, it, mu) = &
         &     dSdpi(ix, iy, it, mu) - tzi_real*sum(aimag(gamval(mu, idirac)* &
                  &(conjg(R(:, ix, iy, it, idirac))* expp * &
                  &        X2(:, ix + ixup, iy + iyup, it + itup, igork1) &
                  &+ conjg(R(:, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                  &             X2(:, ix, iy, it, igork1))))
              enddo
            enddo
          enddo
        endif
        !
      enddo
    enddo

    sumdS_S=sum(dSdpi)
    absdS_S=sum(abs(dSdpi))
    call MPI_AllReduce(MPI_In_Place, sumdS_S, 1, MPI_Real, MPI_Sum, comm, ierr)
    call MPI_AllReduce(MPI_In_Place, absdS_S, 1, MPI_Real, MPI_Sum, comm, ierr)

#ifdef DERIVSWILSON
!      call derivsWilsonDiagonal
!      call derivsWilsonMass
!      call derivsMassTest(R, X2, anum, iflag, am, imass)
!      call derivsDiagTest(R, X2, anum, iflag, am, imass)
      if (iflag.eq.1) then
      call derivsDPTest(R, X2, anum, iflag, am, imass,1,kthird-1)
      call derivsDMTest(R, X2, anum, iflag, am, imass,kthird-1,1)
      call derivsDPTest(R, X2, anum, iflag, am, imass,2,kthird)
      call derivsDMTest(R, X2, anum, iflag, am, imass,kthird,2)
      elseif (iflag.eq.0) then
      call derivsPDTest(R, X2, anum, iflag, am, imass,1,kthird-1)
      call derivsMDTest(R, X2, anum, iflag, am, imass,kthird-1,1)
      call derivsPDTest(R, X2, anum, iflag, am, imass,2,kthird)
      call derivsMDTest(R, X2, anum, iflag, am, imass,kthird,2)
      endif
     sumdS_W=sum(dSdpi)
     absdS_W=sum(abs(dSdpi))
     call MPI_AllReduce(MPI_In_Place, sumdS_W, 1, MPI_Real, MPI_Sum, comm, ierr)
     call MPI_AllReduce(MPI_In_Place, absdS_W, 1, MPI_Real, MPI_Sum, comm, ierr)
#endif

 if (ip_global .eq. 0) then
     open(unit=105,file="fort.105",position="append")
     write(105,'(5f20.12)') sumdS_S,sumdS_W,absdS_S,absdS_W,absdS_W/absdS_S
     close(105)
 endif

      return
      end subroutine derivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsWilsonDiagonal(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass

    real :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, mu
    integer :: igork1
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer idxL,idxR
    integer ierr

    tzi_real = 2*real(anum)
!   print *,"WILSON KERNEL"
    if (COMPACT) then
      print *,"COMPACT NOT IMPLEMENTED FOR WILSON KERNEL"
      stop
    endif

      expp=1d0
      expm=1d0

      ! upper diagonal
      if (iflag.eq.0) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        do idirac=1,4 ! P-
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=1,kthird-1
                  idxL=il
                  idxR=il+1
                  sliceR=cmplx(0.0,0.0)
                  sliceR(:,:,:,3:4)=X2(idxR,:,:,:, 3:4)

                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &     tzi_real*real(akappa)*( aimag( &
                &         conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
                &           sliceR(ix + ixup, iy + iyup, it + itup, idirac)) &
                & - aimag(conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &         X2(idxR, ix, iy, it, idirac)) ) 
                &         sliceR(ix, iy, it, idirac)) ) 
                enddo

              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=1,kthird-1
                  idxL=il
                  idxR=il+1
                  sliceR=cmplx(0.0,0.0)
                  sliceR(:,:,:,3:4)=X2(idxR,:,:,:,3:4)

                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &        tzi_real* aimag( gamval(mu, idirac)* &
                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
                &        sliceR(ix + ixup, iy + iyup, it + itup, igork1) &
                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &             X2(idxR, ix, iy, it, igork1)) )
                &             sliceR(ix, iy, it, igork1)) )
                enddo

              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      elseif (iflag.eq.1) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        do idirac=1,4 ! P+
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=1,kthird-1
                  idxL=il
                  idxR=il+1
                  sliceL=cmplx(0.0,0.0)
                  sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)

                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &     tzi_real*real(akappa)*( aimag( &
!                &         conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &         conjg(sliceL(ix, iy, it, idirac))* expp * &
                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
!                & - aimag(conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                & - aimag(conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &         X2(idxR, ix, iy, it, idirac)) ) 
                enddo

              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=1,kthird-1
                  idxL=il
                  idxR=il+1
                  sliceL=cmplx(0.0,0.0)
                  sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)

                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) - &
                &        tzi_real* aimag( gamval(mu, idirac)* &
!                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &   (conjg(sliceL(ix, iy, it, idirac))* expp * &
                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
!                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &  + conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &             X2(idxR, ix, iy, it, igork1)) )
                enddo

              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      endif ! iflag

      ! lower diagonal
      if (iflag.eq.0) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        do idirac=1,4 ! P+
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=2,kthird
                  idxL=il
                  idxR=il-1
                  sliceR=cmplx(0.0,0.0)
                  sliceR(:,:,:,1:2)=X2(idxR,:,:,:,1:2)

                  dSdpi(ix, iy, it, mu) = &
                &       dSdpi(ix, iy, it, mu) + tzi_real*real(akappa)*( aimag( &
                &         conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
                &           sliceR(ix + ixup, iy + iyup, it + itup, idirac)) &
                & - aimag(conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &         X2(idxR, ix, iy, it, idirac)) ) 
                &         sliceR(ix, iy, it, idirac)) ) 
                enddo
              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=2,kthird
                  idxL=il
                  idxR=il-1
                  sliceR=cmplx(0.0,0.0)
                  sliceR(:,:,:,1:2)=X2(idxR,:,:,:,1:2)

                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &          tzi_real* aimag( gamval(mu, idirac)* &
                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
                &        sliceR(ix + ixup, iy + iyup, it + itup, igork1) &
                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &             X2(idxR, ix, iy, it, igork1)) )
                &             sliceR(ix, iy, it, igork1)) )
                enddo

              enddo
            enddo
          enddo
        end do ! idirac
      end do ! mu

    elseif (iflag.eq.1) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        do idirac=1,4 ! P-
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=2,kthird
                  idxL=il
                  idxR=il-1
                  sliceL=cmplx(0.0,0.0)
                  sliceL(:,:,:,3:4)=R(idxL,:,:,:,3:4)

                  dSdpi(ix, iy, it, mu) = &
                &       dSdpi(ix, iy, it, mu) + tzi_real*real(akappa)*( aimag( &
!                &         conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &         conjg(sliceL(ix, iy, it, idirac))* expp * &
                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
!                & - aimag(conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                & - aimag(conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &         X2(idxR, ix, iy, it, idirac)) ) 
                enddo
              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l

                do il=2,kthird
                  idxL=il
                  idxR=il-1
                  sliceL=cmplx(0.0,0.0)
                  sliceL(:,:,:,3:4)=R(idxL,:,:,:,3:4)

                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) - &
                &          tzi_real* aimag( gamval(mu, idirac)* &
!                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &   (conjg(sliceL(ix, iy, it, idirac))* expp * &
                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
!                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &  + conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &             X2(idxR, ix, iy, it, igork1)) )
                enddo

              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      endif

  return
  end subroutine derivsWilsonDiagonal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsWilsonMass(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass

    real :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, mu
    integer :: igork1
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer idxL,idxR
    real sumdS_W,absdS_W
    integer ierr

    tzi_real = 2*real(anum)

      ! upper right mass 
      if (imass.eq.1) then
        cmult=cmplx(-am,0)
      elseif (imass.eq.3) then
        cmult=cmplx(0,-am)
      endif

      if (iflag.eq.0) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        idxL=1
        idxR=kthird
        sliceR=cmplx(0.0,0.0)
        sliceR(:,:,:,1:2)=X2(idxR,:,:,:,1:2)

        do idirac=1,4 ! P+
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &                 tzi_real*real(akappa)*( aimag( &
                &     cmult* conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
                &           sliceR(ix + ixup, iy + iyup, it + itup, idirac)) &
                & - aimag(cmult* conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &           X2(idxR, ix, iy, it, idirac)) )
                &           sliceR(ix, iy, it, idirac)) )
              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &           tzi_real*aimag( cmult*gamval(mu, idirac)* &
                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
                &        sliceR(ix + ixup, iy + iyup, it + itup, igork1) &
                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &             X2(idxR, ix, iy, it, igork1)) )
                &             sliceR(ix, iy, it, igork1)) )
              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      elseif (iflag.eq.1) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        idxL=1
        idxR=kthird
        sliceL=cmplx(0.0,0.0)
        sliceL(:,:,:,3:4)=R(idxL,:,:,:,3:4)

        do idirac=1,4 ! P-
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &              cmult*tzi_real*real(akappa)*( aimag( &
!                &         conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &       cmult* conjg(sliceL(ix, iy, it, idirac))* expp * &
                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
!                & - aimag(conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                & - aimag(cmult*conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &           X2(idxR, ix, iy, it, idirac)) )
              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) - &
                &            tzi_real*aimag(cmult*gamval(mu, idirac)* &
!                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &   (conjg(sliceL(ix, iy, it, idirac))* expp * &
                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
!                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &  + conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &             X2(idxR, ix, iy, it, igork1)) )
              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      end if

      ! lower left mass 
      if (imass.eq.1) then
        cmult=cmplx(-am,0)
      elseif (imass.eq.3) then
        cmult=cmplx(0,am)
      endif

      if (iflag.eq.0) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        idxL=kthird
        idxR=1
        sliceR=cmplx(0.0,0.0)
        sliceR(:,:,:,3:4)=X2(idxR,:,:,:,3:4)

        do idirac=1,4 ! P-
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &            tzi_real*real(akappa)*( aimag( &
                &       cmult*  conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
                &           sliceR(ix + ixup, iy + iyup, it + itup, idirac)) &
                & - aimag(cmult*conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &           X2(idxR, ix, iy, it, idirac)) )
                &           sliceR(ix, iy, it, idirac)) )
              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &          tzi_real*aimag(cmult*gamval(mu, idirac)* &
                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
!                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
                &        sliceR(ix + ixup, iy + iyup, it + itup, igork1) &
                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
!                &             X2(idxR, ix, iy, it, igork1)) )
                &             sliceR(ix, iy, it, igork1)) )
              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      elseif (iflag.eq.1) then

      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        idxL=kthird
        idxR=1
        sliceL=cmplx(0.0,0.0)
        sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)

        do idirac=1,4 ! P+
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) + &
                &            tzi_real*real(akappa)*( aimag( &
!                &         conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &       cmult*  conjg(sliceL(ix, iy, it, idirac))* expp * &
                &           X2(idxR, ix + ixup, iy + iyup, it + itup, idirac)) &
!                & - aimag(conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                & - aimag(cmult*conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &           X2(idxR, ix, iy, it, idirac)) )
              enddo
            enddo
          enddo

          igork1 = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                  dSdpi(ix, iy, it, mu) = dSdpi(ix, iy, it, mu) - &
                &          tzi_real*aimag(cmult*gamval(mu, idirac)* &
!                &   (conjg(R(idxL, ix, iy, it, idirac))* expp * &
                &   (conjg(sliceL(ix, iy, it, idirac))* expp * &
                &        X2(idxR, ix + ixup, iy + iyup, it + itup, igork1) &
!                &  + conjg(R(idxL, ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &  + conjg(sliceL(ix + ixup, iy + iyup, it + itup, idirac))* expm * &
                &             X2(idxR, ix, iy, it, igork1)) )
              enddo
            enddo
          enddo
        enddo ! idirac
      enddo ! mu

      endif
 
      return
      end subroutine derivsWilsonMass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsDiagTest(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer idxL,idxR

      cmult=cmplx(1.0,0.0)

      ! upper diagonal
!      do il=1,kthird-1
      do il=1,1
        idxL=il
        idxR=il+1
!        sliceL(:,:,:,:)=R(idxL,:,:,:,:)
!        sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        sliceL=cmplx(0.0,0.0)
        sliceR=cmplx(0.0,0.0)
        if (iflag.eq.1) then
          sliceL(:,:,:,:)=R(idxL,:,:,:,:)
          sliceR(:,:,:,1:2)=X2(idxR,:,:,:,3:4)
        elseif (iflag.eq.0) then
          sliceL(:,:,:,3:4)=R(idxL,:,:,:,1:2)
          sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        endif
        call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)
      enddo

!      do il=2,kthird
      do il=2,2
        ! lower left mass 
        idxL=il
        idxR=il-1
        sliceL=cmplx(0.0,0.0)
        sliceR=cmplx(0.0,0.0)
!        sliceL(:,:,:,:)=R(idxL,:,:,:,:)
!        sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        if (iflag.eq.1) then
          sliceL(:,:,:,:)=R(idxL,:,:,:,:)
          sliceR(:,:,:,3:4)=X2(idxR,:,:,:,1:2)
        elseif (iflag.eq.0) then
          sliceL(:,:,:,1:2)=R(idxL,:,:,:,3:4)
          sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        endif
        call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)
      end do

      return
      end subroutine derivsDiagTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsMassTest(R, X2, anum, iflag, am, imass)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer, intent(in) :: iflag,imass

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer idxL,idxR

      cmult=cmplx(-am,0)

        ! upper right mass 
        idxL=1
        idxR=kthird
        sliceL=0.0
        sliceR=0.0
        if (iflag.eq.1) then
          sliceL(:,:,:,:)=R(idxL,:,:,:,:)
          sliceR(:,:,:,1:2)=X2(idxR,:,:,:,1:2)
!          sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        elseif (iflag.eq.0) then
          sliceL(:,:,:,3:4)=R(idxL,:,:,:,3:4)
!          sliceL(:,:,:,:)=R(idxL,:,:,:,:)
          sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        endif
        call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

        ! lower left mass 
        idxL=kthird
        idxR=1
        sliceL=0.0
        sliceR=0.0
        if (iflag.eq.1) then
          sliceL(:,:,:,:)=R(idxL,:,:,:,:)
          sliceR(:,:,:,3:4)=X2(idxR,:,:,:,3:4)
!          sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        elseif (iflag.eq.0) then
          sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)
!          sliceL(:,:,:,:)=R(idxL,:,:,:,:)
          sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
        endif
        call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

      return
      end subroutine derivsMassTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsDTest(R, X2, anum, iflag, am, imass,idxL,idxR)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

      cmult=cmplx(1.0,0.0)
      sliceL(:,:,:,:)=R(idxL,:,:,:,:)
      sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
      call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

      return
      end subroutine derivsDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsDPTest(R, X2, anum, iflag, am, imass,idxL,idxR)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

      cmult=cmplx(1.0,0.0)
      sliceL(:,:,:,:)=R(idxL,:,:,:,:)
      sliceR=cmplx(0.0,0.0)
      sliceR(:,:,:,1:2)=X2(idxR,:,:,:,1:2)
      call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

      return
      end subroutine derivsDPTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsDMTest(R, X2, anum, iflag, am, imass,idxL,idxR)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

      cmult=cmplx(1.0,0.0)
      sliceL(:,:,:,:)=R(idxL,:,:,:,:)
      sliceR=cmplx(0.0,0.0)
      sliceR(:,:,:,3:4)=X2(idxR,:,:,:,3:4)
      call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

      return
      end subroutine derivsDMTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsPDTest(R, X2, anum, iflag, am, imass,idxL,idxR)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

      cmult=cmplx(1.0,0.0)
      sliceL=cmplx(0.0,0.0)
      sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)
      sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
      call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

      return
      end subroutine derivsPDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivsMDTest(R, X2, anum, iflag, am, imass,idxL,idxR)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
    use comms
    implicit none
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: anum
    real,intent(in) :: am
    integer,intent(in) :: iflag,imass
    integer,intent(in) :: idxL,idxR

    integer :: mu
    complex(dp) :: expp,expm,cmult
    real(dp) :: mult
    integer il
    complex(dp) :: sliceL(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: sliceR(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

      cmult=cmplx(1.0,0.0)
      sliceL=cmplx(0.0,0.0)
      sliceL(:,:,:,1:2)=R(idxL,:,:,:,1:2)
      sliceR(:,:,:,:)=X2(idxR,:,:,:,:)
      call derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)

      return
      end subroutine derivsMDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine derivSlice(sliceL,sliceR,anum,iflag,am,imass,cmult)
    use gforce, only: dSdpi
    use gammamatrices, only: kdelta, gamval, gamin
    use params, only: kthird, ksizet_l, ksizey_l, ksizex_l, akappa, dp
    use trial, only: theta
    use params, only: COMPACT
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
    complex(dp) :: expp,expm
    real :: dmult
    integer il
    integer ierr

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

