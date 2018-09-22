module dirac_split
! Versions of the dirac operator subroutines which operate only on a subset 
! of the lattice. 
  use params
  use comms
  implicit none
  save

  complex(dp) :: gamval(6,4)
  integer :: gamin(6,4)

contains 

  pure subroutine dslash_split_nonlocal(Phi,R,u,am,imass,chunk,mu,v)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer :: chunk(2,3) ! portion of array to operate on
    integer, intent(in) :: mu
    integer, intent(in) :: v
    integer :: xd,xu,yd,yu,td,tu ! portion of array to operate on
    integer :: ixup, iyup, itup, ix, iy, it, idirac, igork

    xd=chunk(1,1)
    xu=chunk(2,1)
    yd=chunk(1,2)
    yu=chunk(2,2)
    td=chunk(1,3)
    tu=chunk(2,3)

    ixup = kdelta(1, mu)
    iyup = kdelta(2, mu)
    itup = kdelta(3, mu)

    if(v.eq.1) then
      do idirac=1,4
        igork=gamin(mu,idirac)
        do it = td,tu
          do iy = yd,yu
            do ix = xd,xu
              Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac) &
                ! Wilson term (hermitian)
              &    -akappa*(u(ix,iy,it,mu) &
                &              * R(:, ix+ixup, iy+iyup, it+itup, idirac)) &
                ! Dirac term (antihermitian)
              &     + gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * R(:, ix+ixup, iy+iyup, it+itup, igork) )
            enddo
          enddo
        enddo
      enddo
    else if(v.eq.-1) then
      do idirac=1,4
        igork=gamin(mu,idirac)
        do it = td,tu
          do iy = yd,yu
            do ix = xd,xu
              Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac) &
                ! Wilson term (hermitian)
              &    -akappa*( conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * R(:, ix-ixup, iy-iyup, it-itup, idirac)) &
                ! Dirac term (antihermitian)
              &     + gamval(mu,idirac) * &
                &       (- conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine dslash_split_nonlocal 

  pure subroutine dslash_split_local(am,Phi,R,imass,chunk)
    implicit none
    real, intent(in) :: am
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    integer :: chunk(2,3) ! portion of array to operate on
    integer :: xd,xu,yd,yu,td,tu ! portion of array to operate on
    real :: diag
    complex(dp) :: zkappa

    xd=chunk(1,1)
    xu=chunk(2,1)
    yd=chunk(1,2)
    yu=chunk(2,2)
    td=chunk(1,3)
    tu=chunk(2,3)

    diag=(3.0-am3)+1.0
    Phi(:,xd:xu,yd:yu,td:tu,:) = diag * R(:,xd:xu,yd:yu,td:tu,:)
    !      
    !  s-like term exploiting projection
    Phi(1:kthird-1, xd:xu, yd:yu, td:tu, 3:4) &
      & = Phi(1:kthird-1, xd:xu, yd:yu, td:tu, 3:4) &
      & - R(2:kthird, xd:xu, yd:yu, td:tu, 3:4)
    Phi(2:kthird, xd:xu, yd:yu, td:tu, 1:2) &
      & = Phi(2:kthird, xd:xu, yd:yu, td:tu, 1:2) &
      & - R(1:kthird-1, xd:xu, yd:yu, td:tu, 1:2)
    !
    !  Mass term (couples the two walls unless imass=5)
    if (imass.eq.1) then
      zkappa=cmplx(am,0.0)
      Phi(kthird, xd:xu, yd:yu, td:tu, 3:4) = &
        & Phi(kthird, xd:xu, yd:yu, td:tu, 3:4) &
        & + zkappa * R(1, xd:xu, yd:yu, td:tu, 3:4)
      Phi(1, xd:xu, yd:yu, td:tu, 1:2) = &
        & Phi(1, xd:xu, yd:yu, td:tu, 1:2) + &
        & zkappa * R(kthird, xd:xu, yd:yu, td:tu, 1:2)
    elseif (imass.eq.3) then
      zkappa=cmplx(0.0,-am)
      Phi(kthird, xd:xu, yd:yu, td:tu, 3:4) = &
        & Phi(kthird, xd:xu, yd:yu, td:tu, 3:4) &
        & - zkappa * R(1, xd:xu, yd:yu, td:tu, 3:4)
      Phi(1, xd:xu, yd:yu, td:tu, 1:2) = &
        & Phi(1, xd:xu, yd:yu, td:tu, 1:2) &
        & + zkappa * R(kthird, xd:xu, yd:yu, td:tu, 1:2)
    elseif (imass.eq.5) then
      zkappa=cmplx(0.0,-am)
      Phi(kthird, xd:xu, yd:yu, td:tu, 3:4) = &
        & Phi(kthird, xd:xu, yd:yu, td:tu, 3:4) &
        & - zkappa * R(kthird, xd:xu, yd:yu, td:tu, 1:2)
      Phi(1, xd:xu, yd:yu, td:tu, 1:2) = &
        & Phi(1, xd:xu, yd:yu, td:tu, 1:2) &
        & - zkappa * R(1, xd:xu, yd:yu, td:tu, 3:4)
    endif
    !
    return
  end subroutine dslash_split_local

  pure subroutine dslashd_split_nonlocal(Phi,R,u,am,imass,chunk,mu,v)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on
    integer, intent(in) :: mu
    integer, intent(in) :: v
    integer :: xd,xu,yd,yu,td,tu ! portion of array to operate on
    integer :: ixup, iyup, itup, ix, iy, it, idirac, igork

    xd=chunk(1,1)
    xu=chunk(2,1)
    yd=chunk(1,2)
    yu=chunk(2,2)
    td=chunk(1,3)
    tu=chunk(2,3)

    ixup = kdelta(1, mu)
    iyup = kdelta(2, mu)
    itup = kdelta(3, mu)

    if(v.eq.1) then
      do idirac=1,4
        igork=gamin(mu,idirac)
        do it = td,tu 
          do iy = yd,yu 
            do ix = xd,xu
              phi(:,ix,iy,it,idirac)=phi(:,ix,iy,it,idirac) &
                !   wilson term (hermitian)
              &    - akappa * (u(ix,iy,it,mu) &
                &              * r(:, ix+ixup, iy+iyup, it+itup, idirac)) &
                !   dirac term (antihermitian)
              &    - gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * r(:, ix+ixup, iy+iyup, it+itup, igork))
            enddo
          enddo
        enddo
      enddo
    else if(v.eq.-1) then
      do idirac=1,4
        igork=gamin(mu,idirac)
        do it = td,tu 
          do iy = yd,yu 
            do ix = xd,xu
              phi(:,ix,iy,it,idirac)=phi(:,ix,iy,it,idirac) &
                !   wilson term (hermitian)
              &    - akappa * (conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * r(:, ix-ixup, iy-iyup, it-itup, idirac)) &
                !   dirac term (antihermitian)
              &    - gamval(mu,idirac) * &
                &       (- conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * r(:, ix-ixup, iy-iyup, it-itup, igork))
            enddo
          enddo
        enddo
      enddo
    endif


  end subroutine dslashd_split_nonlocal

  pure subroutine dslashd_split_local(am,Phi,R,imass,chunk)
    implicit none
    real, intent(in) :: am
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    integer :: chunk(2,3) ! portion of array to operate on
    integer :: xd,xu,yd,yu,td,tu ! portion of array to operate on
    real :: diag
    complex(dp) :: zkappa

    xd=chunk(1,1)
    xu=chunk(2,1)
    yd=chunk(1,2)
    yu=chunk(2,2)
    td=chunk(1,3)
    tu=chunk(2,3)

    diag=(3.0-am3)+1.0
    Phi(:,xd:xu,yd:yu,td:tu,:) = diag * R(:,xd:xu,yd:yu,td:tu,:)

    !   s-like term exploiting projection
    Phi(1:kthird-1, xd:xu, yd:yu, td:tu, 1:2) &
      & = Phi(1:kthird-1, xd:xu, yd:yu, td:tu, 1:2) &
      & - R(2:kthird, xd:xu, yd:yu, td:tu, 1:2)
    Phi(2:kthird, xd:xu, yd:yu, td:tu, 3:4) &
      & = Phi(2:kthird, xd:xu, yd:yu, td:tu, 3:4) &
      & - R(1:kthird-1, xd:xu, yd:yu, td:tu, 3:4)
    !
    !   Mass term (couples the two walls unless imass=5) 
    if(imass.eq.1)then
      zkappa=cmplx(am,0.0)
      Phi(kthird, xd:xu, yd:yu, td:tu, 1:2) = &
        & Phi(kthird, xd:xu, yd:yu, td:tu, 1:2) &
        & + zkappa * R(1, xd:xu, yd:yu, td:tu, 1:2)
      Phi(1, xd:xu, yd:yu, td:tu, 3:4) = &
        & Phi(1, xd:xu, yd:yu, td:tu, 3:4) &
        & + zkappa * R(kthird, xd:xu, yd:yu, td:tu, 3:4)
    elseif(imass.eq.3)then
      zkappa = cmplx(0.0,am)
      Phi(kthird, xd:xu, yd:yu, td:tu, 1:2) = &
        & Phi(kthird, xd:xu, yd:yu, td:tu, 1:2) &
        & + zkappa * R(1, xd:xu, yd:yu, td:tu, 1:2)
      Phi(1, xd:xu, yd:yu, td:tu, 3:4) = &
        & Phi(1, xd:xu, yd:yu, td:tu, 3:4) &
        & - zkappa * R(kthird, xd:xu, yd:yu, td:tu, 3:4)
    elseif(imass.eq.5)then
      zkappa = cmplx(0.0,am)
      Phi(kthird, xd:xu, yd:yu, td:tu, 1:2) = &
        & Phi(kthird, xd:xu, yd:yu, td:tu, 1:2) &
        & - zkappa * R(kthird, xd:xu, yd:yu, td:tu, 3:4)
      Phi(1, xd:xu, yd:yu, td:tu, 3:4) = &
        & Phi(1,xd:xu, yd:yu, td:tu, 3:4) &
        & - zkappa * R(1, xd:xu, yd:yu, td:tu, 1:2)
    endif

    return 
  end subroutine dslashd_split_local

  pure integer function kdelta(nu, mu)
    integer, intent(in) :: nu
    integer, intent(in) :: mu

    kdelta=merge(1,0,nu==mu)
  end function kdelta

end module dirac_split





