module diracWilson
  use params
  use gammamatrices
  implicit none
  save
  real(dp) :: coeffs(kthird)
  real(dp) :: revcoeffs(kthird)

contains
 
      pure subroutine DWilson(Phi,R,u,mass)
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
      pure subroutine DWilsonD(Phi,R,u,mass)
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
  pure subroutine dslash(Phi,R,u,am,imass)
    implicit none
    !     calculates Phi = M*R
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer il

    Phi=cmplx(0.0,0.0)

    ! P+ component
    Rslice=cmplx(0.0,0.0)
    do il=1,kthird
      Rslice(:,:,:,1:2)=R(il,:,:,:,1:2)
      call DWilson(Phislice,Rslice,u,-am3)
      ! diagonal
      Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+Phislice+Rslice
      ! lower diagonal
      if (il.lt.kthird) then
        Phi(il+1,:,:,:,:)=Phi(il+1,:,:,:,:)+Phislice-Rslice
      end if
    end do

    ! P- component
    Rslice=cmplx(0.0,0.0)
    do il=1,kthird
      Rslice(:,:,:,3:4)=R(il,:,:,:,3:4)
      call DWilson(Phislice,Rslice,u,-am3)
      ! diagonal
      Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+Phislice+Rslice
      ! upper diagonal
      if (il.gt.1) then
        Phi(il-1,:,:,:,:)=Phi(il-1,:,:,:,:)+Phislice-Rslice
      end if
    end do

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
  end subroutine dslash
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure subroutine dslashslow(Phi,R,u,am,imass)
    implicit none
    !     calculates Phi = M*R
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!    real(dp),intent(in) :: coeffs(kthird)
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
      Phi(il,:,:,:,:)=coeffs(il)*Phislice+Rslice
    end do

      !  upper diagonal
      do il=1,kthird-1
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,3:4)=R(il+1,:,:,:,3:4)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    + coeffs(il)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do

      ! lower diagonal
      do il=2,kthird
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,1:2)=R(il-1,:,:,:,1:2)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &  +   coeffs(il)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do
      
      !  Mass term - top right
      Rslice=cmplx(0.0,0.0)
      Rslice(:,:,:,1:2)=R(kthird,:,:,:,1:2)
      call DWilson(Phislice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:) = -am*(coeffs(1)*Phislice-Rslice)
      elseif (imass.eq.3) then
        Mslice(:,:,:,:) = -cmplx(0.0,am)*(coeffs(1)*Phislice-Rslice)
      endif
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &                Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  & +             Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

      !  Mass term - bottom left
      Rslice=cmplx(0.0,0.0)
      Rslice(:,:,:,3:4)=R(1,:,:,:,3:4)
      call DWilson(Phislice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:)=-am*(coeffs(kthird)*Phislice-Rslice)
      elseif (imass.eq.3) then
        Mslice(:,:,:,:)=cmplx(0.0,am)*(coeffs(kthird)*Phislice-Rslice)
      endif
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &                 Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            +        Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)


    !
    return
  end subroutine dslashslow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dslashd(Phi,R,u,am,imass,reqs_R)
    use comms, only : complete_halo_update
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer, dimension(12),intent(inout), optional :: reqs_R
    integer il
    complex(dp) :: zkappa
    real :: diag

    Phi=cmplx(0.0,0.0)
    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+Rslice
      ! lower diagonal
      if (il.lt.kthird) then
        Phi(il+1,:,:,:,3:4) = Phi(il+1,:,:,:,3:4)   &
 &    - Rslice(:,:,:,3:4)
      end if
      ! upper diagonal
      if (il.gt.1) then
        Phi(il-1,:,:,:,1:2) = Phi(il-1,:,:,:,1:2)   &
 &    - Rslice(:,:,:,1:2)
      end if
    enddo


    if(present(reqs_r)) then
      call complete_halo_update(reqs_R)
    endif

    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+Phislice
      ! lower diagonal
      if (il.lt.kthird) then
        Phi(il+1,:,:,:,3:4) = Phi(il+1,:,:,:,3:4)   &
 &    + Phislice(:,:,:,3:4)
      end if
      ! upper diagonal
      if (il.gt.1) then
        Phi(il-1,:,:,:,1:2) = Phi(il-1,:,:,:,1:2)   &
 &    + Phislice(:,:,:,1:2)
      end if
    enddo

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
  end subroutine dslashd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dslashdslow(Phi,R,u,am,imass,reqs_R)
    use comms, only : complete_halo_update
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
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
      Phi(il,:,:,:,:)=revcoeffs(il)*Phislice+Rslice
    enddo

      ! upper diagonal
      do il=1,kthird-1
        Rslice=R(il+1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &  +   revcoeffs(il+1)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)
      enddo

      ! lower diagonal
      do il=2,kthird
        Rslice=R(il-1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)   &
 &  +   revcoeffs(il-1)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)
      enddo

      ! upper right mass term  
      Rslice=R(kthird,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Mslice=revcoeffs(kthird)*Phislice-Rslice
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
      Mslice=revcoeffs(1)*Phislice-Rslice
      if (imass.eq.1) then
        zkappa = -cmplx(am,0.0)
      elseif (imass.eq.3) then
        zkappa = cmplx(0.0,am)
      endif
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
  &                 Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) &
  &            + zkappa*Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)

    return
  end subroutine dslashdslow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dslashShamirMass(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag

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
  end subroutine dslashShamirMass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dslashdShamirMass(Phi,R,u,am,imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag

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
      print *,"imass=5 for dslashdShamirMass no implemented"
      stop
    endif

    return
  end subroutine dslashdShamirMass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepZolo()
  use zolomodule
  use ratfuncs
  implicit none

  real(dp) lmin,lmax
  type(zolotarev) :: zolo
  type(sgnratfunc) :: srf
  integer i

  lmin=1d-3
  lmax=20
  call setZolo(lmin,lmax,kthird,zolo)
  call getRoots(zolo) 
  coeffs=1d0/zolo%roots
!  print *,"coeffs",coeffs
  do i=1,kthird
    revcoeffs(i)=coeffs(kthird-i+1)
!    revcoeffs(i)=coeffs(i)
  end do
!  call testZolo(28) 

  return
  end subroutine prepZolo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module diracWilson

