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

    ! diagonal
    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(il,:,:,:,:)=Phislice+Rslice
    end do

!    call dslash_diagShamir(Phi,R,am,imass)
!    call dslash_massShamir(Phi,R,am,imass)

    call dslash_diagWilson(Phi,R,u,am,imass)
    call dslash_massWilson(Phi,R,u,am,imass)

!    call dslash_diagTest(Phi,R,u,am,imass)
!    call dslash_massTest(Phi,R,u,am,imass)

!    call dslash_DTest(Phi,R,u,am,imass,1,kthird-1)
!    call dslash_DTest(Phi,R,u,am,imass,kthird-1,1)
!    call dslash_DTest(Phi,R,u,am,imass,2,kthird)
!    call dslash_DTest(Phi,R,u,am,imass,kthird,2)

!    call dslash_DMTest(Phi,R,u,am,imass,1,kthird-1)
!    call dslash_DPTest(Phi,R,u,am,imass,kthird-1,1)
!    call dslash_DMTest(Phi,R,u,am,imass,2,kthird)
!    call dslash_DPTest(Phi,R,u,am,imass,kthird,2)

!    call dslash_DMTest(Phi,R,u,am,imass,1,2)
!    call dslash_DPTest(Phi,R,u,am,imass,2,1)
!    call dslash_DMTest(Phi,R,u,am,imass,kthird-1,kthird)
!    call dslash_DPTest(Phi,R,u,am,imass,kthird,kthird-1)

!    do il=1,kthird-1
!    call dslash_DMTest(Phi,R,u,am,imass,il,il+1)
!    call dslash_DPTest(Phi,R,u,am,imass,il+1,il)
!    enddo

    return
  end subroutine dslash
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_diagShamir(Phi,R,am,imass)
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am

    !  upper diagonal
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

    ! lower diagonal
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

    return
  end subroutine dslash_diagShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_diagWilson(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer il

      !  upper diagonal
      do il=1,kthird-1
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,3:4)=R(il+1,:,:,:,3:4)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do

      ! lower diagonal
      do il=2,kthird
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,1:2)=R(il-1,:,:,:,1:2)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do
      
    return
  end subroutine dslash_diagWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_diagTest(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer il

      !  upper diagonal
!      do il=1,kthird-1
      do il=1,1
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,3:4)=R(il+1,:,:,:,3:4)
!        Rslice(:,:,:,:)=R(il+1,:,:,:,:)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do

      ! lower diagonal
!      do il=2,kthird
      do il=2,2
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,1:2)=R(il-1,:,:,:,1:2)
!        Rslice(:,:,:,:)=R(il-1,:,:,:,:)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do
      
    return
  end subroutine dslash_diagTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_massShamir(Phi,R,am,imass)
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa

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
  end subroutine dslash_massShamir

    subroutine dslash_massWilson(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer il

      !  Mass term - top right
      Rslice=cmplx(0.0,0.0)
      Rslice(:,:,:,1:2)=R(kthird,:,:,:,1:2)
      call DWilson(Phislice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:) = -am*(Phislice-Rslice)
      elseif (imass.eq.3) then
        Mslice(:,:,:,:) = -cmplx(0.0,am)*(Phislice-Rslice)
      endif
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &                Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  & +             Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

      !  Mass term - bottom left
      Rslice=cmplx(0.0,0.0)
      Rslice(:,:,:,3:4)=R(1,:,:,:,3:4)
      call DWilson(Phislice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:)=-am*(Phislice-Rslice)
      elseif (imass.eq.3) then
        Mslice(:,:,:,:)=cmplx(0.0,am)*(Phislice-Rslice)
      endif
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &                 Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            +        Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

    !
    return
  end subroutine dslash_massWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_massTest(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer idx
    complex(dp) :: su(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1)

      ! upper right mass term - 1:2
      Rslice=0
      Rslice(:,:,:,1:2)=R(kthird,:,:,:,1:2)
!      Rslice(:,:,:,:)=R(kthird,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            - am * Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

      ! lower left mass term - 3:4 
      Rslice=0
      Rslice(:,:,:,3:4)=R(1,:,:,:,3:4)
!      Rslice(:,:,:,:)=R(1,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            - am * Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslash_massTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_DTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice(:,:,:,:)=R(idxR,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslash_DTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_DPTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice=0
      Rslice(:,:,:,1:2)=R(idxR,:,:,:,1:2)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslash_DPTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_DMTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice=0
      Rslice(:,:,:,3:4)=R(idxR,:,:,:,3:4)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) & 
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslash_DMTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_PDTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice(:,:,:,:)=R(idxR,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) 

    return
    end subroutine dslash_PDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslash_MDTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice(:,:,:,:)=R(idxR,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) 

    return
    end subroutine dslash_MDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!    call dslashd_diagShamir(Phi,R,am,imass)
!    call dslashd_massShamir(Phi,R,am,imass)

    call dslashd_diagWilson(Phi,R,u,am,imass)
    call dslashd_massWilson(Phi,R,u,am,imass)

!    call dslashd_diagTest(Phi,R,u,am,imass)
!    call dslashd_massTest(Phi,R,u,am,imass)

!    call dslashd_DTest(Phi,R,u,am,imass,1,kthird-1)
!    call dslashd_DTest(Phi,R,u,am,imass,kthird-1,1)
!    call dslashd_DTest(Phi,R,u,am,imass,2,kthird)
!    call dslashd_DTest(Phi,R,u,am,imass,kthird,2)

!    call dslashd_PDTest(Phi,R,u,am,imass,1,kthird-1)
!    call dslashd_MDTest(Phi,R,u,am,imass,kthird-1,1)
!    call dslashd_PDTest(Phi,R,u,am,imass,2,kthird)
!    call dslashd_MDTest(Phi,R,u,am,imass,kthird,2)

!    call dslashd_PDTest(Phi,R,u,am,imass,1,2)
!    call dslashd_MDTest(Phi,R,u,am,imass,2,1)
!    call dslashd_PDTest(Phi,R,u,am,imass,kthird-1,kthird)
!    call dslashd_MDTest(Phi,R,u,am,imass,kthird,kthird-1)

!    do il=1,kthird-1
!    call dslashd_PDTest(Phi,R,u,am,imass,il,il+1)
!    call dslashd_MDTest(Phi,R,u,am,imass,il+1,il)
!    enddo

    return
  end subroutine dslashd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dslashd_diagShamir(Phi,R,am,imass)
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am

    !   s-like term exploiting projection
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

    return
  end subroutine dslashd_diagShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_diagWilson(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer il

      ! upper diagonal
      do il=1,kthird-1
        Rslice=R(il+1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)
      enddo

      ! lower diagonal
      do il=2,kthird
        Rslice=R(il-1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)
      enddo

    return
  end subroutine dslashd_diagWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_diagTest(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer il

      ! upper diagonal - 1:2
!      do il=1,kthird-1
      do il=1,1
        Rslice=R(il+1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)
      enddo

      ! lower diagonal - 3:4
!      do il=2,kthird
      do il=2,2
        Rslice=R(il-1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)   &
 &  + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)
      enddo

    return
  end subroutine dslashd_diagTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dslashd_massShamir(Phi,R,am,imass)
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa

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
    endif

    return
  end subroutine dslashd_massShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_massWilson(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer il
    complex(dp) :: zkappa

      ! upper right mass term  
      Rslice=R(kthird,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Mslice=Phislice-Rslice
      if (imass.eq.1) then
        zkappa = -cmplx(am,0.0)
      elseif (imass.eq.3) then
        zkappa = -cmplx(0.0,am)
      endif
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
  &                Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
  & +      zkappa*Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)

      ! lower left mass term  
      Rslice=R(1,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Mslice=Phislice-Rslice
      if (imass.eq.1) then
        zkappa = -cmplx(am,0.0)
      elseif (imass.eq.3) then
        zkappa = cmplx(0.0,am)
      endif
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
  &                 Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) &
  &            + zkappa*Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)

    return
  end subroutine dslashd_massWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_massTest(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer idx
    complex(dp) :: su(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1)

      ! upper right mass term - 3:4 
      Rslice=R(kthird,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
  &   Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) &
  &            - am * Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) 

      ! lower left mass term - 1:2
      Rslice=R(1,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
  &   Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) &
  &            - am * Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) 

    return
    end subroutine dslashd_massTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_DTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice=R(idxR,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslashd_DTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_DPTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice=0
      Rslice(:,:,:,1:2)=R(idxR,:,:,:,1:2)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslashd_DPTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_DMTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice=0
      Rslice(:,:,:,3:4)=R(idxR,:,:,:,3:4)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) & 
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:) 

    return
    end subroutine dslashd_DMTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_PDTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice(:,:,:,:)=R(idxR,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) 

    return
    end subroutine dslashd_PDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd_MDTest(Phi,R,u,am,imass,idxL,idxR)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer,intent(in) :: idxL,idxR
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      Rslice(:,:,:,:)=R(idxR,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
  &   Phi(idxL,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
  &            + Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) &
  &            - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) 

    return
    end subroutine dslashd_MDTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module diracShamir

