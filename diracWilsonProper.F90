module diracWilson
  use params
  use gammamatrices
  implicit none
  save
  contains
 
    pure subroutine DWilson(Phi,R,u,mass)
    implicit none
    complex(dp), intent(out) :: Phi(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: mass
    real diag
    integer mu,idirac,it,iy,ix,ixup,iyup,itup,igork

    ! diagonal term
    diag=(3.0+mass)! 3 is from Wilson term
    Phi=diag*R

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

      ! diagonal
      do il=1,kthird
        Rslice=R(il,:,:,:,:)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,:,:,:,:)=Phislice+Rslice
      end do

      return

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

      return
      end subroutine dslash
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dslashd(Phi,R,u,am,imass,reqs_R)
    use comms, only : complete_halo_update
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

      ! diagonal
      do il=1,kthird
        Rslice=R(il,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,:,:,:,:)=Phislice+Rslice
      enddo

      return

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine diagV(Xin,Xout,u)
      complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
      complex(dp), intent(out) :: Xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(in) :: Xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      integer il

      ! diagonal
      do il=1,kthird
        Rslice=Xin(il,:,:,:,:)
        call DWilson(Phislice,Rslice,u,-am3)
        Xout(il,:,:,:,:)=Phislice-Rslice
      end do

      return
      end subroutine diagV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reflect5(Xin,Xout)
      implicit none
      complex(dp), intent(in) :: Xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(out) :: Xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      integer il

      ! diagonal
      do il=1,kthird
        Xout(il,:,:,:,:)=Xin(kthird+1-il,:,:,:,:)
      end do

      return
      end subroutine reflect5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mGam4(Xin,Xout,gid)
      use gammamatrices
      use params
      implicit none
      complex(dp) :: Xin(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: Xout(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      integer gid,igork
      integer idirac,it,iy,ix
      complex(dp) val

      do idirac=1,4
        igork=gamin(gid,idirac)
        val=gamval(gid,idirac)/akappa
        do it = 1,ksizet_l
          do iy = 1,ksizey_l
            do ix = 1,ksizex_l
              Xout(ix,iy,it,idirac)=val*Xin(ix,iy,it,igork) 
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine mGam4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mGam5(Xin,Xout,gid)
      use gammamatrices
      use params
      implicit none
      complex(dp) :: Xin(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: Xout(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      integer gid,igork
      integer idirac,it,iy,ix,il
      complex(dp) val

      do idirac=1,4
        igork=gamin(gid,idirac)
        val=gamval(gid,idirac)/akappa
        do it = 1,ksizet_l
          do iy = 1,ksizey_l
            do ix = 1,ksizex_l
              do il=1,kthird
                Xout(il,ix,iy,it,idirac)=val*Xin(il,ix,iy,it,igork) 
              enddo
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine mGam5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine test
      use dum1
      use gaussian
      use comms
      use comms4
      use comms5
      implicit none
      complex(dp) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
      complex(dp) :: R4(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: R4init(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: TMP4(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: D4(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: G3D4G3(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: R5(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: R5init(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: TMP5(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: D5(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: G3D5G3(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: LHS(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: RHS(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      integer :: reqs_ps(12),reqs_TMP5(12),reqs_TMP4(12),reqs_R4(12),reqs_R5(12)
      integer, dimension(12) :: reqs_u
      real mass
      integer it,iy,ix,idirac,mu,il

      mass=0.1   
   
!     setup u
      do mu=1,3
      call gauss0(ps, reqs_ps)
      call complete_halo_update(reqs_ps)
      do it = 1,ksizet_l
        do iy = 1,ksizey_l
          do ix = 1,ksizex_l
            u(ix,iy,it,mu) = 1.0 + cmplx(0.0,ps(ix,iy,it,1))
          end do
        end do
      end do
      end do
      if (ibound .eq. -1 .and. ip_t .eq. (np_t - 1)) then
        u(:, :, ksizet_l, 3) = -u(:, :, ksizet_l, 3)
      end if
      call start_halo_update_4(3, u, 3, reqs_u)
      call complete_halo_update(reqs_u)

!     setup Rinit
      do il=1,kthird
      do idirac=1,4
      call gauss0(ps, reqs_ps)
      call complete_halo_update(reqs_ps)
      do it = 1,ksizet_l
        do iy = 1,ksizey_l
          do ix = 1,ksizex_l
            R4(ix,iy,it,idirac)=cmplx(ps(ix,iy,it,1),ps(ix,iy,it,2))
            R5(il,ix,iy,it,idirac)=cmplx(ps(ix,iy,it,1),ps(ix,iy,it,2))
          end do
        end do
      end do
      end do
      end do
      call start_halo_update_4(4, R4, 2, reqs_R4)
      call complete_halo_update(reqs_R4)
      R4init=R4
      call start_halo_update_5(4, R5, 3, reqs_R5)
      call complete_halo_update(reqs_R5)
      R5init=R5

      ! calc G3.D.G3.R
      call mGam4(R4,TMP4,4)
      call mGam5(R5,TMP5,4)

!      print *,"TMP4:",TMP4(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

      call start_halo_update_4(4, TMP4, 4, reqs_TMP4)
      call complete_halo_update(reqs_TMP4)
      call DWilson(R4,TMP4,u,mass)

      call start_halo_update_5(4, TMP5, 5, reqs_TMP5)
      call complete_halo_update(reqs_TMP5)
      call dslash(R5,TMP5,u,mass,1)

      call start_halo_update_4(4, R4, 6, reqs_R4)
      call complete_halo_update(reqs_R4)
      call  mGam4(R4,G3D4G3,4)

      call start_halo_update_5(4, R5, 7, reqs_R5)
      call complete_halo_update(reqs_R5)
      call  mGam5(R5,G3D5G3,4)


      ! calc Ddag.R
      call  DWilsonD(D4,R4init,u,mass)
      call dslashd(D5,R5init,u,mass,1)

      print *,"max err:",maxval(abs(D4-G3D4G3))
      print *,"max err:",maxval(abs(D4(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)-G3D4G3(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)))

      print *,"max err:",maxval(abs(D5-G3D5G3))
      print *,"max err:",maxval(abs(D5(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)-G3D5G3(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)))
      
!     dslash.V.R.g5.R = V.R.g5.dshashd.R

      call mGam5(R5init,TMP5,4)
      call reflect5(TMP5,R5)
      call start_halo_update_5(4, R5, 11, reqs_R5)
      call complete_halo_update(reqs_R5)
      call diagV(R5,TMP5,u)
      call start_halo_update_5(4, TMP5, 12, reqs_TMP5)
      call complete_halo_update(reqs_TMP5)
      call dslash(LHS,TMP5,u,mass,1)


      call dslashd(TMP5,R5init,u,mass,1)
      call start_halo_update_5(4, TMP5, 13, reqs_TMP5)
      call complete_halo_update(reqs_TMP5)
      call mGam5(TMP5,R5,4)
      call reflect5(R5,TMP5)
      call start_halo_update_5(4, TMP5, 14, reqs_TMP5)
      call complete_halo_update(reqs_TMP5)
      call diagV(TMP5,RHS,u)

      print *,"max err:",maxval(abs(LHS-RHS))
      print *,"max err:",maxval(abs(LHS(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)-RHS(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)))

      return
      end subroutine test

end module diracWilson

