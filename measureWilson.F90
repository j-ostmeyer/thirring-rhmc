module measureWilson
  implicit none
  contains
  !******************************************************************
  !    matrix inversion via conjugate gradient algorithm
  !       solves (Mdagger)Mx=Phi,
  !           NB. no even/odd partitioning
  !******************************************************************
  subroutine congrad(Phi, res, itercg, am, imass, iterations)
    use trial, only: u
    use vector
    use comms5, only: init_halo_update_5,start_halo_update_5
    use comms_common, only: comm
    use comms
    use diracWilson
    implicit none
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    integer, intent(in) :: imass
    integer, intent(out), optional :: iterations
    complex(dp) :: x1(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: x2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: p(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: r(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real :: resid
    real(dp) :: betacg, betacgn, betacgd, alpha, alphan, alphad
    integer :: i
    integer, dimension(12) :: reqs_x1, reqs_r
    integer :: ierr

    resid = 4*ksize*ksize*ksizet*kthird*res*res
    itercg = 0
    alphan = 0.0

!    call init_halo_update_5(4, x1, 8, reqs_x1) ! ATTEMPT
!    call init_halo_update_5(4, r, 10, reqs_r)  ! ATTEMPT

!   initialise
    r = Phi ! = RR is input
    x=r ! (DR=RR) is output

    call dslash(x1, r, u, am, imass) ! call Mptr(RR,x1,u,DAGGER,mass,add)
!    call MPI_Startall(12, reqs_x1, ierr)
    call start_halo_update_5(4, x1, 8, reqs_x1) ! ATTEMPT
    call complete_halo_update(reqs_x1)
    call dslashd(x2, x1, u, am, imass) ! call Mptr(x1,x2,u,.not.DAGGER,mass,add)
    r=Phi-x2 ! r=RR-x2
!    call MPI_Startall(12, reqs_r, ierr)
    call start_halo_update_5(4, r, 10, reqs_r) ! ATTEMPT
    call complete_halo_update(reqs_r)
    p=r
    betacgn=sum(conjg(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
    call MPI_AllReduce(MPI_In_Place, betacgn, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
    if (betacgn.lt.resid) goto 8
      
    do i=1,niterc
        itercg=itercg+1
        call dslash(x1, p, u, am, imass) ! call Mptr(p,x1,u,DAGGER,mass,add)
!        call MPI_Startall(12, reqs_x1, ierr)
        call start_halo_update_5(4, x1, 28, reqs_x1) ! ATTEMPT
!        call complete_halo_update(reqs_x1)
        call dslashd(x2, x1, u, am, imass,reqs_x1) ! call Mptr(x1,x2,u,.not.DAGGER,mass,add)
        alphan=sum(conjg(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
        call MPI_AllReduce(MPI_In_Place, alphan, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
        alphad=sum(conjg(p(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*x2(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
        call MPI_AllReduce(MPI_In_Place, alphad, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
        alpha=alphan/alphad

        x=x+alpha*p
        r=r-alpha*x2
!        call MPI_Startall(12, reqs_r, ierr)
        call start_halo_update_5(4, r, 30, reqs_r) ! ATTEMPT
!        call complete_halo_update(reqs_r)
        betacgn=sum(conjg(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
        call MPI_AllReduce(MPI_In_Place, betacgn, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
        call complete_halo_update(reqs_r)

        if (betacgn.lt.resid) goto 8
        betacg=betacgn/alphan
        p=r+betacg*p
    end do

8     continue

    if (present(iterations)) then
      iterations = i
    endif

    return
  end subroutine congrad

  !*****************************************************************
  !   Calculate fermion expectation values via a noisy estimator
  !   -matrix inversion via conjugate gradient algorithm
  !       solves Mx=x1
  !     (Numerical Recipes section 2.10 pp.70-73)
  !*******************************************************************
  subroutine measureW(psibarpsi, res, aviter, am, imass, isweep_total)
    use trial, only: u
    use vector, xi => x
    use comms5, only: start_halo_update_5
    use comms
    use gaussian
    use params
    use diracWilson
    implicit none
    real, intent(out) :: psibarpsi, aviter
    real, intent(in) :: res, am
    integer, intent(in) :: imass
    integer, intent(in), optional :: isweep_total
    integer, parameter :: knoise = 100
    complex(dp) :: x(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: psibarpsi1, psibarpsi2
    complex(dp) :: oslice(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: islice(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp) :: cnum(0:1), cden(1)
    real :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
    real :: pt(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
    real(dp) :: pbp(knoise)
    integer :: idsource, idsource2, inoise
    integer :: iter, itercg
    real :: susclsing
    integer, dimension(12) :: reqs_ps, reqs_pt, reqs_Phi
    integer :: ierr

    print *,"Wilson measure"

    iter = 0
    cnum(0) = 0.0
    cnum(1) = 1.0
    cden(1) = 0.0

    do inoise = 1, knoise
!      print *,"inoise:",inoise
      !     set up noise

      call gauss0(ps, reqs_ps)
!DB      call complete_halo_update(reqs_ps)
!ORIG      psibarpsi1 = (0.0, 0.0)
      call gauss0(pt, reqs_pt)
!DB      call complete_halo_update(reqs_pt)
!ORIG      psibarpsi2 = (0.0, 0.0)
      
        psibarpsi1 = (0.0, 0.0)
        psibarpsi2 = (0.0, 0.0)
!ORIG      do idsource = 1, 2
      do idsource = 1, 4
!        call gauss0(ps, reqs_ps)
!        call gauss0(pt, reqs_pt)
        !
        !  source on domain wall at ithird=1
        !
        x = cmplx(0.0, 0.0)
        !  We started a halo update in the gauss0 call;
        !  Now we need it to be complete if it isn't already
        if (idsource .eq. 1) then
          call complete_halo_update(reqs_ps)
        end if

        x(:, :, :, idsource) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))
        !
        xi = cmplx(0.0, 0.0)
        xi(1, :, :, :, :) = x
        !
        ! Phi= Mdagger*xi
        !
        call dslashd(Phi, xi, u, am, imass)
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, Phi, 16, reqs_Phi)
        call complete_halo_update(reqs_Phi)

        call congrad(Phi, res, itercg, am, imass)
        iter = iter + itercg
        !
! ORIG        ! for idirac=1 or 2 
! JW    ! for P+
        islice=cmplx(0,0);
        islice(:,:,:,1:2)=xi(kthird,:,:,:,1:2) ! P+.phi
        call DWilson(oslice,islice,u,-am3)
        oslice=oslice-islice
        xi(kthird,:,:,:,:)=oslice
          
        psibarpsi1 = psibarpsi1 &
    &           + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))* &
    &               xi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))
      end do
        !  source on domain wall at ithird=kthird
!ORIG        idsource2 = idsource + 2
!        idsource2 = idsource + 2
        idsource2 = idsource 
      do idsource2 = 1,4
!        call gauss0(pt, reqs_pt)
        x = cmplx(0.0, 0.0)
        !  Again, if this isn't finished by now we have to wait for it
!ORIG        if (idsource .eq. 1) then
        if (idsource2 .eq. 1) then
          call complete_halo_update(reqs_pt)
        end if

        x(:, :, :, idsource2) = cmplx(pt(:, :, :, 1), pt(:, :, :, 2))
        xi = cmplx(0.0, 0.0)
        xi(kthird, :, :, :, :) = x
        call dslashd(Phi, xi, u, am, imass)

        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, Phi, 18, reqs_Phi)
        call complete_halo_update(reqs_Phi)

        ! xi= (M)**-1 * Phi
        !
        call congrad(Phi, res, itercg, am, imass)
        iter = iter + itercg
        !
        ! for idirac=3 or 4
        islice=cmplx(0,0);
        islice(:,:,:,3:4)=xi(1,:,:,:,3:4) ! P-.phi
        call DWilson(oslice,islice,u,-am3)
        oslice=oslice-islice
        xi(1,:,:,:,:)=oslice

          ! pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
          psibarpsi2 = psibarpsi2 &
    &           + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2)) &
        &                *xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2))
        !  end trace on Dirac indices....
      enddo ! idsource

      !  The sums psibarpsi1 and psibarpsi2 are initialised to 0 outside the loop
      !  So can be summed up per process, then collected here at the end
      call MPI_AllReduce(MPI_In_Place, psibarpsi1, 1, MPI_Double_Complex, &
        & MPI_Sum, comm, ierr)
      call MPI_AllReduce(MPI_In_Place, psibarpsi2, 1, MPI_Double_Complex, &
        & MPI_Sum, comm, ierr)

      !
      if (imass .eq. 1) then
        psibarpsi1 = psibarpsi1/kvol
        psibarpsi2 = psibarpsi2/kvol
        pbp(inoise) = psibarpsi1 + psibarpsi2
      elseif (imass .eq. 3) then
        psibarpsi1 = cmplx(0.0, +1.0)*psibarpsi1/kvol
        psibarpsi2 = cmplx(0.0, -1.0)*psibarpsi2/kvol
        pbp(inoise) = psibarpsi1 + psibarpsi2
      endif

      if (ip_global .eq. 0) then
        open (unit=100, file='fort.100', action='write', position='append')
        if (present(isweep_total)) then
          write (100, '(I5,4E19.9E3)') isweep_total, real(psibarpsi1), aimag(psibarpsi1), &
             & real(psibarpsi2), aimag(psibarpsi2)
        else
          write (100, '(4E19.9E3)') real(psibarpsi1), aimag(psibarpsi1), &
             & real(psibarpsi2), aimag(psibarpsi2)
        endif
        close (100)
      end if

      ! end loop on noise
    enddo
!      print *,ip_global,"loop ended"
    !
    susclsing = 0.0
    !
    psibarpsi = real(sum(pbp))
!      print *,ip_global,"loop ended"
    do inoise = 1, knoise
      susclsing = susclsing + real(sum(pbp(inoise)*pbp(inoise + 1:knoise)))
    enddo
    psibarpsi = psibarpsi/knoise
    susclsing = 2*kvol*susclsing/(knoise*(knoise - 1))
    if (ip_global .eq. 0) then
!      print *,"print fort.200"
      open (unit=200, file='fort.200', action='write', position='append')
      if (present(isweep_total)) then
        write (200, '(I5,3E17.7E3)') isweep_total, psibarpsi, susclsing
      else
        write (200, '(2E17.7E3)') psibarpsi, susclsing
      endif
      close (200)
    end if
    aviter = float(iter)/(4*knoise)
    return
  end subroutine measureW

end module measureWilson
