module diracShamir
  use params
  use gammamatrices
  implicit none
  save

contains
 
      subroutine DWilson(Phi,R,u,mass)
      implicit none
      complex(dp), intent(out) :: Phi(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
      real, intent(in) :: mass
      real(dp) diag
      integer mu,idirac,it,iy,ix,ixup,iyup,itup,igork

    print *,"DWilson in diracShamir"
    stop

    ! diagonal term
    diag=(3.0+mass)! 3 is from Wilson term
    Phi=diag*R
    !      
    !     Wilson term (hermitian) and Dirac term (antihermitian)
    do mu=1,3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4
        igork=gamin(mu,idirac)
        do it = 1,ksizet_l
          do iy = 1,ksizey_l
            do ix = 1,ksizex_l
              Phi(ix,iy,it,idirac)=Phi(ix,iy,it,idirac) &
                ! Wilson term (hermitian)
              &    -akappa*(u(ix,iy,it,mu) &
                &              * R(ix+ixup, iy+iyup, it+itup, idirac) &
                &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * R(ix-ixup, iy-iyup, it-itup, idirac)) &
               ! Dirac term (antihermitian)
              &     + gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * R(ix+ixup, iy+iyup, it+itup, igork) &
                &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * R(ix-ixup, iy-iyup, it-itup, igork))
            enddo
          enddo
        enddo
      enddo
    enddo

      return
      end subroutine DWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWilsonD(Phi,R,u,mass)
      implicit none
      complex(dp), intent(out) :: Phi(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
      real, intent(in) :: mass
      integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
      real diag

      diag=(3.0+mass)
      Phi=diag*R
      do mu=1,3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        do idirac=1,4
          igork=gamin(mu,idirac)
          do it = 1,ksizet_l
            do iy = 1,ksizey_l
              do ix = 1,ksizex_l
                phi(ix,iy,it,idirac)=phi(ix,iy,it,idirac) &
                !   wilson term (hermitian)
                &    - akappa * (u(ix,iy,it,mu) &
                &              * r(ix+ixup, iy+iyup, it+itup, idirac) &
                &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * r(ix-ixup, iy-iyup, it-itup, idirac)) &
                !   dirac term (antihermitian)
                &    - gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * r(ix+ixup, iy+iyup, it+itup, igork) &
                &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * r(ix-ixup, iy-iyup, it-itup, igork))
              enddo
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine DWilsonD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dslash(Phi,R,u,am,imass)
    !     calculates Phi = M*R
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer il

    print *,"new Dirac Shamir"
    stop

    ! diagonal
    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(il,:,:,:,:)=Phislice+Rslice
    end do

    !  upper diagonal
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

    ! lower diagonal
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

    !  Mass term 
    if (imass.eq.1) then
      zkappa=cmplx(am,0.0)
      Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        & + zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) + &
        & zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    elseif (imass.eq.3) then
      zkappa=cmplx(0.0,-am)
      Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        & - zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        & + zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    endif
    !
    return
  end subroutine dslash
  !***********************************************************************
  subroutine dslashd(Phi,R,u,am,imass,reqs_R)
    use comms, only : complete_halo_update
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer, dimension(12),intent(inout), optional :: reqs_R
    integer il
    complex(dp) :: zkappa
    real :: diag

    if(present(reqs_r)) then
      call complete_halo_update(reqs_R)
    endif

    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(il,:,:,:,:)=Phislice+Rslice
    enddo

    !   s-like term exploiting projection
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    !
    !   Mass term (couples the two walls unless imass=5) 
    if(imass.eq.1)then
      zkappa=cmplx(am,0.0)
      Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        & + zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        & + zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    elseif(imass.eq.3)then
      zkappa = cmplx(0.0,am)
      Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        & + zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        & - zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    elseif(imass.eq.5)then
    endif

    return
  end subroutine dslashd

end module diracShamir

