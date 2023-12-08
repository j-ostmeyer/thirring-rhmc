module dirac
  use comms_common
  use params
  use gammamatrices

  implicit none
  save

contains

  pure subroutine dslash(Phi, R, u, am, imass)
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am

    call verify_kernel_choice()

#ifdef SHAMIR_KERNEL
    call dslash_shamir(Phi, R, u, am, imass)
#endif

#ifdef WILSON_KERNEL
    call dslash_wilson(Phi, R, u, am, imass)
#endif

  end subroutine

  pure subroutine dslash_shamir(Phi, R, u, am, imass)
    !
    !     calculates Phi = M*R
    !
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(out) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    !     write(6,*) 'hi from dslash'
    !
    !     diagonal term
    diag = (3.0 - am3) + 1.0
    Phi = diag*R
    !
    !     Wilson term (hermitian) and Dirac term (antihermitian)
    do mu = 1, 3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac = 1, 4
        igork = gamin(mu, idirac)
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              Phi(1:kthird_l, ix, iy, it, idirac) = Phi(1:kthird_l, ix, iy, it, idirac) &
                                                    ! Wilson term (hermitian)
                                                    - akappa &
                                                    *(u(ix, iy, it, mu) &
                                                      *R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac) &
                                                      + conjg(u(ix - ixup, iy - iyup, it - itup, mu)) &
                                                      *R(1:kthird_l, ix - ixup, iy - iyup, it - itup, idirac)) &
                                                    ! Dirac term (antihermitian)
                                                    + gamval(mu, idirac) &
                                                    *(u(ix, iy, it, mu) &
                                                      *R(1:kthird_l, ix + ixup, iy + iyup, it + itup, igork) &
                                                      - conjg(u(ix - ixup, iy - iyup, it - itup, mu)) &
                                                      *R(1:kthird_l, ix - ixup, iy - iyup, it - itup, igork))
            enddo
          enddo
        enddo
      enddo
    enddo

    ! s-like term exploiting projection
    ! If only one rank along third dimension
    if (np_third .eq. 1) then
      Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        - R(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

      Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        - R(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      ! If more than one rank along third dimension
    else
      if (ip_third .eq. 0) then
        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - R(2:kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

        Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - R(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

      else if (ip_third .eq. np_third - 1) then
        Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - R(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - R(0:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

      else
        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - R(2:kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - R(0:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      end if
    end if

    ! Mass term (couples the two walls unless imass=5)
    if (imass .eq. 1) then
      zkappa = cmplx(am, 0.0)

      ! Separate if statements. Needed when np_third == 1
      if (ip_third .eq. np_third - 1) then
        Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          + zkappa*R(kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      end if
      if (ip_third .eq. 0) then
        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          + zkappa*R(0, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      end if

    else if (imass .eq. 3) then
      zkappa = cmplx(0.0, -am)

      ! Separate if statements. Needed when np_third == 1
      if (ip_third .eq. np_third - 1) then
        Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - zkappa*R(kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      end if
      if (ip_third .eq. 0) then
        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          + zkappa*R(0, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      end if

    else if (imass .eq. 5) then
      zkappa = cmplx(0.0, -am)
      !         do idirac=3,4
      !         igork=gamin(5,idirac)

      ! Separate if statements. Needed when np_third == 1
      if (ip_third .eq. np_third - 1) then
        Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - zkappa*R(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
        !        Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) = &
        !            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) &
        !            & + 2 * zkappa * gamval(5,idirac) * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, igork)
        !         enddo
        !         do idirac=1,2
        !         igork=gamin(5,idirac)
      end if
      if (ip_third .eq. 0) then
        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - zkappa*R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
        !        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) = &
        !            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac)
        !            & + 2 * zkappa * gamval(5,idirac) * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, igork)
        !         enddo
      end if
    end if

    return
  end subroutine dslash_shamir

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

  pure subroutine dslash_wilson(Phi,R,u,am,imass)
    implicit none
    !     calculates Phi = M*R
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(0:kthird+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:kthird+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
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
  end subroutine dslash_wilson

  !***********************************************************************

  pure subroutine dslashd_local_shamir(am, Phi, R, imass)
    implicit none
    complex(dp), intent(out) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real :: diag
    complex(dp) :: zkappa
    real, intent(in) :: am
    integer, intent(in) :: imass

    diag = (3.0 - am3) + 1.0
    Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = diag*R(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)

    ! s-like term exploiting projection
    ! If only one rank along third dimension
    if (np_third == 1) then
      Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
        Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
        - R(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

      Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
        Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
        - R(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      ! If more than one rank along third dimension
    else
      if (ip_third .eq. 0) then
        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - R(2:kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

        Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - R(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

      else if (ip_third .eq. np_third - 1) then
        Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - R(2:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - R(0:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)

      else
        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - R(2:kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)

        Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - R(0:kthird_l - 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      end if
    end if

    ! Mass term (couples the two walls unless imass=5)
    if (imass .eq. 1) then
      zkappa = cmplx(am, 0.0)

      ! Separate if statements. Needed when np_third == 1
      if (ip_third .eq. np_third - 1) then
        Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          + zkappa*R(kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      end if
      if (ip_third .eq. 0) then
        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          + zkappa*R(0, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      end if

    else if (imass .eq. 3) then
      zkappa = cmplx(0.0, am)

      ! Separate if statements. Needed when np_third == 1
      if (ip_third .eq. np_third - 1) then
        Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          + zkappa*R(kthird_l + 1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      end if
      if (ip_third .eq. 0) then
        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - zkappa*R(0, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      end if

    else if (imass .eq. 5) then
      zkappa = cmplx(0.0, am)

      ! Separate if statements. Needed when np_third == 1
      if (ip_third .eq. np_third - 1) then
        Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
          Phi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
          - zkappa*R(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
      end if
      if (ip_third .eq. 0) then
        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
          Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
          - zkappa*R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
      end if

    end if

    return
  end subroutine dslashd_local_shamir

!***********************************************************************

#ifdef MPI
  subroutine dslashd(Phi, R, u, am, imass, reqs_R)
#else
  pure subroutine dslashd(Phi, R, u, am, imass, reqs_R)
#endif
    use comms, only: complete_halo_update
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer, dimension(16), intent(inout), optional :: reqs_R

    call verify_kernel_choice()

#ifdef SHAMIR_KERNEL
  call dslashd_shamir(Phi, R, u, am, imass, reqs_R)
#endif

#ifdef WILSON_KERNEL
  call dslashd_wilson(Phi, R, u, am, imass, reqs_R)
#endif

end subroutine dslashd

#ifdef MPI
  subroutine dslashd_shamir(Phi, R, u, am, imass, reqs_R)
#else
  pure subroutine dslashd_shamir(Phi, R, u, am, imass, reqs_R)
#endif
      use comms, only: complete_halo_update
      implicit none
      complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
      complex(dp), intent(out) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
      complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
      integer, intent(in) :: imass
      real, intent(in) :: am
      integer, dimension(16), intent(inout), optional :: reqs_R

      integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork

      ! We need to update the halo before calling dslashd_local
      ! because the halo is now needed along the "third" dimension
#ifdef MPI
      if (present(reqs_R)) then
        call complete_halo_update(reqs_R)
      end if
#endif

      !   taking care of the part that does not need the halo (not valid anymore-halo is needed)
      !   diagonal term (hermitian)
      call dslashd_local_shamir(am, Phi, R, imass)
      !   call complete_halo_update_5(4, Phi)

      !   taking care of the part that does need the halo
      !   wilson term (hermitian) and dirac term (antihermitian)
      do mu = 1, 3
        ixup = kdelta(1, mu)
        iyup = kdelta(2, mu)
        itup = kdelta(3, mu)

        do idirac = 1, 4
          igork = gamin(mu, idirac)
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                Phi(1:kthird_l, ix, iy, it, idirac) = Phi(1:kthird_l, ix, iy, it, idirac) &
                                                      ! Wilson term (hermitian)
                                                      - akappa &
                                                      *(u(ix, iy, it, mu) &
                                                        *R(1:kthird_l, ix + ixup, iy + iyup, it + itup, idirac) &
                                                        + conjg(u(ix - ixup, iy - iyup, it - itup, mu)) &
                                                        *R(1:kthird_l, ix - ixup, iy - iyup, it - itup, idirac)) &
                                                      ! Dirac term (antihermitian)
                                                      - gamval(mu, idirac) &
                                                      *(u(ix, iy, it, mu) &
                                                        *R(1:kthird_l, ix + ixup, iy + iyup, it + itup, igork) &
                                                        - conjg(u(ix - ixup, iy - iyup, it - itup, mu)) &
                                                        *R(1:kthird_l, ix - ixup, iy - iyup, it - itup, igork))
              enddo
            enddo
          enddo
        enddo
      enddo

      return
  end subroutine dslashd_shamir

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

#ifdef MPI
  subroutine dslashd_wilson(Phi, R, u, am, imass, reqs_R)
#else
  pure subroutine dslashd_wilson(Phi, R, u, am, imass, reqs_R)
#endif
    use comms, only : complete_halo_update
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(0:kthird+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:kthird+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer, dimension(16),intent(inout), optional :: reqs_R
    integer il
    complex(dp) :: zkappa
    real :: diag

#ifdef MPI
      if (present(reqs_R)) then
        call complete_halo_update(reqs_R)
      end if
#endif

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
  end subroutine dslashd_wilson

  !***********************************************************************
  pure subroutine dslash2d(Phi, R, u)
    !     calculates Phi = m*R
    !     complex, intent(in) :: u(ksizex_l, ksizey_l, ksizet_l, 3)
    !     complex, intent(out) :: Phi(ksizex_l,ksizey_l,ksizet_l, 4)
    !     complex, intent(in) :: R(ksizex_l,ksizey_l,ksizet_l,4)
    implicit none
    complex(dp), intent(in) ::  u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(out) :: Phi(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer :: ix, iy, it, idirac, mu, ixup, iyup, igork
    real(dp) :: diag

    !     diagonal term
    diag = 2.0d0
    Phi = diag*R

    !     wilson and dirac terms
    do mu = 1, 2
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)

      do idirac = 1, 4
        igork = gamin(mu, idirac)
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              Phi(ix, iy, it, idirac) = Phi(ix, iy, it, idirac) &
                                        ! wilson term
                                        - akappa &
                                        *(u(ix, iy, it, mu) &
                                          *R(ix + ixup, iy + iyup, it, idirac) &
                                          + conjg(u(ix - ixup, iy - iyup, it, mu)) &
                                          *R(ix - ixup, iy - iyup, it, idirac)) &
                                        ! dirac term
                                        + gamval(mu, idirac) &
                                        *(u(ix, iy, it, mu) &
                                          *R(ix + ixup, iy + iyup, it, igork) &
                                          - conjg(u(ix - ixup, iy - iyup, it, mu)) &
                                          *R(ix - ixup, iy - iyup, it, igork))
            enddo
          enddo
        enddo
      enddo
    enddo
    !    call complete_halo_update_4(4, Phi)

    return
  end subroutine dslash2d
  !***********************************************************************
  !   A Kronecker delta function
  !   Useful for calculating coordinate offsets
  !***********************************************************************

end module dirac
