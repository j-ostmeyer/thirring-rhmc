module diracWilson
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

      diag=(3.0-mass)
      Phi=diag*R
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
      real(dp) diag

      diag=(3.0-mass)
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
      complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
      integer, intent(in) :: imass
      real, intent(in) :: am
      complex(dp) :: Rmp(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp),dimension(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4) :: Rslice,Mslice
      complex(dp) :: zkappa
      real :: diag
      integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork, il

!     diagonal term
      diag=(3.0-am3)+1.0
      Phi=diag*R
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
                Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac) &
                ! Wilson term (hermitian)
                &    -akappa*(u(ix,iy,it,mu) &
                &              * R(:, ix+ixup, iy+iyup, it+itup, idirac) &
                &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * R(:, ix-ixup, iy-iyup, it-itup, idirac)) &
                ! Dirac term (antihermitian)
                &     + gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * R(:, ix+ixup, iy+iyup, it+itup, igork) &
                &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
              enddo
            enddo
          enddo
        enddo
      enddo

      diag=-diag
      Rmp=0d0
      Rmp(:,:,:,:,3:4)=R(:,:,:,:,3:4)
      ! upper diagonal
      do il=1,kthird-1
        Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+diag*Rmp(il+1,:,:,:,:)
        ! Wilson term (hermitian) and Dirac term (antihermitian)
        do mu=1,3
          ixup = kdelta(1, mu)
          iyup = kdelta(2, mu)
          itup = kdelta(3, mu)

          do idirac=1,4
            igork=gamin(mu,idirac)
            do it = 1,ksizet_l
              do iy = 1,ksizey_l
                do ix = 1,ksizex_l
                  Phi(il,ix,iy,it,idirac)=Phi(il,ix,iy,it,idirac) &
                  ! Wilson term (hermitian)
                  &    -akappa*(u(ix,iy,it,mu) &
                  &              * Rmp(il+1, ix+ixup, iy+iyup, it+itup, idirac) &
                  &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &              * Rmp(il+1, ix-ixup, iy-iyup, it-itup, idirac)) &
                  ! Dirac term (antihermitian)
                  &     + gamval(mu,idirac) * &
                  &       (u(ix,iy,it,mu) &
                  &         * Rmp(il+1, ix+ixup, iy+iyup, it+itup, igork) &
                  &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &         * Rmp(il+1, ix-ixup, iy-iyup, it-itup, igork))
                enddo
              enddo
            enddo
          enddo
        enddo
      end do

      Rmp=0d0
      Rmp(:,:,:,:,1:2)=R(:,:,:,:,1:2)
      ! upper diagonal
      do il=2,kthird
        Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+diag*Rmp(il-1,:,:,:,:)
        ! Wilson term (hermitian) and Dirac term (antihermitian)
        do mu=1,3
          ixup = kdelta(1, mu)
          iyup = kdelta(2, mu)
          itup = kdelta(3, mu)
  
          do idirac=1,4
            igork=gamin(mu,idirac)
            do it = 1,ksizet_l
              do iy = 1,ksizey_l
                do ix = 1,ksizex_l
                  Phi(il,ix,iy,it,idirac)=Phi(il,ix,iy,it,idirac) &
                  ! Wilson term (hermitian)
                  &    -akappa*(u(ix,iy,it,mu) &
                  &              * Rmp(il-1, ix+ixup, iy+iyup, it+itup, idirac) &
                  &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &              * Rmp(il-1, ix-ixup, iy-iyup, it-itup, idirac)) &
                  ! Dirac term (antihermitian)
                  &     + gamval(mu,idirac) * &
                  &       (u(ix,iy,it,mu) &
                  &         * Rmp(il-1, ix+ixup, iy+iyup, it+itup, igork) &
                  &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &         * Rmp(il-1, ix-ixup, iy-iyup, it-itup, igork))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      !  Mass term - top right
      Rslice=0d0
      Mslice=0d0
      Rslice(:,:,:,1:2)=Phi(kthird,:,:,:,1:2)
      call DWilson(Mslice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:)=-am*(Mslice-Rslice)
      elseif (imass.eq.3) then
        Mslice(:,:,:,:)=-cmplx(0,am)*(Mslice-Rslice)
      endif
      Phi(1,:,:,:,:)=Phi(1,:,:,:,:)+Mslice

      !  Mass term - bottom left
      Rslice=0d0
      Mslice=0d0
      Rslice(:,:,:,3:4)=Phi(1,:,:,:,3:4)
      call DWilson(Mslice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:)=-am*(Mslice-Rslice)
      elseif (imass.eq.3) then
        Mslice(:,:,:,:)=cmplx(0,am)*(Mslice-Rslice)
      endif
      Phi(kthird,:,:,:,:)=Phi(kthird,:,:,:,:)+Mslice

      return
      end subroutine dslash
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pure subroutine dslashd_local(am,Phi,R,imass)
      implicit none
      complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      real :: diag
      complex(dp) :: zkappa
      real, intent(in) :: am
      integer, intent(in) :: imass

      diag=(3.0-am3)+1.0
      Phi(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = diag * R(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

      return 
      end subroutine dslashd_local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
      subroutine dslashd(Phi,R,u,am,imass,reqs_R)
#else
      subroutine dslashd(Phi,R,u,am,imass)
#endif
      use comms, only : complete_halo_update
      complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
      complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      integer, intent(in) :: imass
      real, intent(in) :: am
      integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork, il
      real(dp) diag
      complex(dp) :: slice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
#ifdef MPI
      integer, dimension(12),intent(inout), optional :: reqs_R
#endif
      !   taking care of the part that does not need the halo
      !   diagonal term (hermitian)
      !call dslashd_local(am,Phi,R,imass)

      !   taking care of the part that does need the halo
      !   wilson term (hermitian) and dirac term (antihermitian)
#ifdef MPI
      if(present(reqs_r)) then
        call complete_halo_update(reqs_R)
      endif
#endif

      Phi=0d0
      ! 1st slice
      call DWilsonD(slice,R(1,:,:,:,:),u,-am3)
      Phi(1,:,:,:,:)=slice+R(1,:,:,:,:)
      Phi(2,:,:,:,1:2)=slice(:,:,:,1:2)-R(1,:,:,:,1:2)
      
      ! middle slices
      do il=2,kthird-1
        call DWilsonD(slice,R(il,:,:,:,:),u,-am3)
        Phi(il-1,:,:,:,3:4)=Phi(il-1,:,:,:,3:4)+ slice(:,:,:,3:4)-R(il,:,:,:,3:4)
        Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+ slice(:,:,:,:)+R(il,:,:,:,:)
        Phi(il+1,:,:,:,1:2)=Phi(il+1,:,:,:,1:2)+ slice(:,:,:,1:2)-R(il,:,:,:,1:2)
      end do

      ! kthird slice
      il=kthird
      call DWilsonD(slice,R(il,:,:,:,:),u,-am3)
      Phi(il-1,:,:,:,3:4)=Phi(il-1,:,:,:,3:4)+ slice(:,:,:,3:4)-R(il,:,:,:,3:4)
      Phi(il,:,:,:,:)=Phi(il,:,:,:,:)+ slice(:,:,:,:)+R(il,:,:,:,:)


      ! upper right mass
      if (imass.eq.1) then
        Phi(1,:,:,:,1:2)=Phi(1,:,:,:,1:2) - am*(slice(:,:,:,1:2)-R(kthird,:,:,:,1:2))
      elseif (imass.eq.3) then
        Phi(1,:,:,:,1:2)=Phi(1,:,:,:,1:2) - cmplx(0,am)*(slice(:,:,:,1:2)-R(kthird,:,:,:,1:2))
      endif

      ! lower left mass
      call DWilsonD(slice,R(1,:,:,:,:),u,-am3)
      if (imass.eq.1) then
        Phi(kthird,:,:,:,3:4)=Phi(kthird,:,:,:,3:4) - am*(slice(:,:,:,3:4)-R(1,:,:,:,3:4))
      elseif (imass.eq.3) then
        Phi(kthird,:,:,:,3:4)=Phi(kthird,:,:,:,3:4) - cmplx(0,am)*(slice(:,:,:,3:4)-R(2,:,:,:,3:4))
      endif

      return
      end subroutine dslashd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !***********************************************************************
  pure subroutine dslash2d(phi,r,u)
    !     calculates phi = m*r
    !
    !     complex, intent(in) :: u(ksizex_l, ksizey_l, ksizet_l, 3)
    !     complex, intent(out) :: phi(ksizex_l,ksizey_l,ksizet_l, 4)
    !     complex, intent(in) :: r(ksizex_l,ksizey_l,ksizet_l,4)
    complex(dp), intent(in) ::  u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: phi(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: r(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer :: ix, iy, it, idirac, mu, ixup, iyup, igork
    real(dp) :: diag

    !     diagonal term
    diag=2.0d0
    phi = diag * r

    !     wilson and dirac terms
    do mu=1,2
      ixup=kdelta(1,mu)
      iyup=kdelta(2,mu)
      !
      do idirac=1,4
        igork=gamin(mu,idirac)
        do it=1,ksizet_l
          do iy=1,ksizey_l
            do ix=1,ksizex_l
              phi(ix,iy,it,idirac) = &
                ! wilson term
              &    phi(ix,iy,it,idirac) &
                &    - akappa * (u(ix,iy,it,mu) * r(ix+ixup, iy+iyup, it, idirac) &
                &             + conjg(u(ix-ixup, iy-iyup, it, mu)) &
                &              * r(ix-ixup, iy-iyup, it, idirac)) &
                ! dirac term
              &     + gamval(mu,idirac) * &
                &      (u(ix,iy,it,mu)*r(ix+ixup, iy+iyup, it, igork) &
                &       - conjg(u(ix-ixup, iy-iyup, it,mu)) &
                &        * r(ix-ixup, iy-iyup, it, igork))
            enddo
          enddo
        enddo
      enddo
    enddo
    !    call complete_halo_update_4(4, phi)
    !
    return
  end subroutine dslash2d
  !!***********************************************************************
  !   A Kronecker delta function
  !   Useful for calculating coordinate offsets
  !***********************************************************************

end module diracWilson

