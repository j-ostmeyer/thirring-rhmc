module dirac
  use params
  use gammamatrices
  implicit none
  save

contains 

  subroutine dslash(Phi,R,u,am,imass)
    !
    !     calculates Phi = M*R
    !
    !     complex, intent(in) :: u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
    !     complex, intent(in) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
    !     complex, intent(in) :: R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
    !     complex :: zkappa
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    !     write(6,*) 'hi from dslash'
    
!    print *,"original Dirac Shamir"
!    stop

    !     diagonal term
    diag=(3.0-am3)+1.0
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
    !
    !  s-like term exploiting projection
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
      & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
      & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    !
    !  Mass term (couples the two walls unless imass=5)
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
    elseif (imass.eq.5) then
      zkappa=cmplx(0.0,-am)
      !         do idirac=3,4
      !         igork=gamin(5,idirac)
      Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        & - zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      !        Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) = &
      !            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) &
      !            & + 2 * zkappa * gamval(5,idirac) * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, igork)
      !         enddo
      !         do idirac=1,2
      !         igork=gamin(5,idirac)
      Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        & - zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      !        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) = &
      !            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) 
      !            & + 2 * zkappa * gamval(5,idirac) * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, igork)
      !         enddo
    endif
    !
    return
  end subroutine dslash
  !***********************************************************************
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
      zkappa = cmplx(0.0,am)
      Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        & - zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        & Phi(1,1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        & - zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    endif

    return 
  end subroutine dslashd_local
#ifdef MPI
  subroutine dslashd(Phi,R,u,am,imass,reqs_R)
#else
  pure subroutine dslashd(Phi,R,u,am,imass)
#endif
    use comms, only : complete_halo_update
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    !complex(dp) :: zkappa
    !real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
#ifdef MPI
    integer, dimension(12),intent(inout), optional :: reqs_R
#endif
    !
    !   taking care of the part that does not need the halo
    !   diagonal term (hermitian)
    call dslashd_local(am,Phi,R,imass)
    !   call complete_halo_update_5(4, phi)
    !
    !   taking care of the part that does need the halo
    !   wilson term (hermitian) and dirac term (antihermitian)
#ifdef MPI
    if(present(reqs_r)) then
      call complete_halo_update(reqs_R)
    endif
#endif
    do mu=1,3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4
        igork=gamin(mu,idirac)
        do it = 1,ksizet_l
          do iy = 1,ksizey_l
            do ix = 1,ksizex_l
              phi(:,ix,iy,it,idirac)=phi(:,ix,iy,it,idirac) &
                !   wilson term (hermitian)
              &    - akappa * (u(ix,iy,it,mu) &
                &              * r(:, ix+ixup, iy+iyup, it+itup, idirac) &
                &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * r(:, ix-ixup, iy-iyup, it-itup, idirac)) &
                !   dirac term (antihermitian)
              &    - gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * r(:, ix+ixup, iy+iyup, it+itup, igork) &
                &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * r(:, ix-ixup, iy-iyup, it-itup, igork))
            enddo
          enddo
        enddo
      enddo
    enddo
    !
    return
#ifdef MPI
  end subroutine dslashd
#else
  end subroutine dslashd
#endif


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

end module dirac

