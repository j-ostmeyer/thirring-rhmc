module measure_module

contains
  !******************************************************************
  !    matrix inversion via conjugate gradient algorithm
  !       solves (Mdagger)Mx=Phi,
  !           NB. no even/odd partitioning
  !******************************************************************
  subroutine congrad(Phi, res, itercg, am, imass, iterations)
    use trial, only: u
    use vector
    use comms5, only: init_halo_update_5
    use comms_common, only: comm
    use comms
    use dirac
    use params
    implicit none
    complex(dp), intent(in) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    integer, intent(in) :: imass
    integer, intent(out), optional :: iterations

    !     complex x1(kferm),x2(kferm),p(kferm),r(kferm)
    complex(dp) :: x1(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: x2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: p(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: r(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real :: resid
    real(dp) :: betacg, betacgn, betacgd, alpha, alphan, alphad
    integer :: nx
#ifdef MPI
    integer, dimension(16) :: reqs_x1, reqs_r
    integer :: ierr
!#else
    !    integer :: reqs_x1, reqs_r
#endif
    !    write(6,111)
    !111 format(' Hi from congrad')
    !
    resid = 4*ksize*ksize*ksizet*kthird*res*res
    itercg = 0
    alphan = 0.0
    !   initialise p=x, r=Phi(na)
    p = x
    r = Phi
    !
    betacgd = 1.0
    alpha = 1.0
    !
#ifdef MPI
    call init_halo_update_5(4, x1, 8, reqs_x1) ! ATTEMPT
    call init_halo_update_5(4, r, 10, reqs_r)  ! ATTEMPT
#endif

    do nx = 1, niterc
      itercg = itercg + 1
      !
      !  x1=Mp
      call dslash(x1, p, u, am, imass)
#ifdef MPI
      call MPI_Startall(16, reqs_x1, ierr)
      !call start_halo_update_5(4, x1, 8, reqs_x1) ! ATTEMPT
#endif
      !
      if (nx .ne. 1) then
        !
        !   alpha=(r,r)/(p,(Mdagger)Mp)
        !   Don't need x1's halo at this point
        alphad = sum(abs(x1(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
#ifdef MPI
        call MPI_AllReduce(MPI_In_Place, alphad, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
#endif
        alpha = alphan/alphad
        !
        !   x=x+alpha*p
        x = x + alpha*p
      end if

      !   Now we need x1's halo, so ensure communication is finished
#ifdef MPI
      call complete_halo_update(reqs_x1)
#else
      call update_halo_5(4, x1)
#endif
      !
      !   x2=(Mdagger)x1=(Mdagger)Mp
      call dslashd(x2, x1, u, am, imass)
      !
      !   r=r-alpha*(Mdagger)Mp
      !   Use x2 with wrong halo since we can't hide communications here
      r = r - alpha*x2

      !   Now update halo for r instead since we can hide communication during the summation
      !   x2 is discarded so we no longer care about its halo
#ifdef MPI
      call MPI_Startall(16, reqs_r, ierr)              ! ATTEMPT
      !call start_halo_update_5(4, r, 10, reqs_r) ! ATTEMPT
#endif

      !   betacg=(r_k+1,r_k+1)/(r_k,r_k)
      betacgn = sum(abs(r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
#ifdef MPI
      call MPI_AllReduce(MPI_In_Place, betacgn, 1, MPI_Double_precision, MPI_Sum, comm, ierr)
#endif
      betacg = betacgn/betacgd
      betacgd = betacgn
      alphan = betacgn

      if (nx .eq. 1) betacg = 0.0

      !   p=r+betacg*p
      !   Now the correct value of r is needed to avoid having to communicate p as well
#ifdef MPI
      call complete_halo_update(reqs_r) ! ATTEMPT
#else
      call update_halo_5(4, r)
#endif
      p = r + betacg*p
      if (betacgn .lt. resid) exit
    end do
    !     write(6,1000)

    if (nx .gt. niterc) then
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
        write (7, 1000)
1000    format(' # iterations of congrad exceeds niterc')
        write (7, *) "Iterations:", nx

#ifdef MPI
      end if
#endif
    endif
    if (present(iterations)) then
      iterations = nx
    endif
    return
  end subroutine congrad

  subroutine measure(psibarpsi, res, aviter, am, imass, isweep_total)
    implicit none
    real, intent(out) :: psibarpsi, aviter
    real, intent(in) :: res, am
    integer, intent(in) :: imass
    integer, intent(in), optional :: isweep_total
    
#ifdef MEASURE_SHAMIR
    call measure_shamir(psibarpsi, res, aviter, am, imass, isweep_total)
#endif

#ifdef MEASURE_WILSON
    call measure_wilson(psibarpsi, res, aviter, am, imass, isweep_total)
#endif
  end subroutine measure

  !*****************************************************************
  !   Calculate fermion expectation values via a noisy estimator
  !   -matrix inversion via conjugate gradient algorithm
  !       solves Mx=x1
  !     (Numerical Recipes section 2.10 pp.70-73)
  !*******************************************************************
  subroutine measure_shamir(psibarpsi, res, aviter, am, imass, isweep_total)
    use trial, only: u
    use vector, xi => x
    use comms5, only: start_halo_update_5
    use comms
    use gaussian
    use dirac
    !use qmrherm_module
    real, intent(out) :: psibarpsi, aviter
    real, intent(in) :: res, am
    integer, intent(in) :: imass
    integer, intent(in), optional :: isweep_total
    !     complex :: x(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    !     complex :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    !     complex :: psibarpsi1,psibarpsi2
    complex(dp) :: x(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: psibarpsi1, psibarpsi2
    real(dp) :: cnum(0:1), cden(1)
    real :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
    real :: pt(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
    real(dp) :: pbp(knoise)
    integer :: idsource, idsource2, inoise
    integer :: iter, itercg
    real :: susclsing
#ifdef MPI
    integer, dimension(16) :: reqs_Phi, reqs_xi
    integer, dimension(12) :: reqs_ps, reqs_pt
    integer :: ierr
#endif
    !     write(6,*) 'hi from measure'
    !
    iter = 0
    !     pbp=0.0
    cnum(0) = 0.0
    cnum(1) = 1.0
    cden(1) = 0.0
    !
    do inoise = 1, knoise
      !
      !     set up noise
#ifdef MPI
      if (ip_third .eq. 0) then
        call gauss0(ps, reqs_ps)
        call gauss0(pt, reqs_pt)
      endif
      psibarpsi1 = (0.0, 0.0)
      psibarpsi2 = (0.0, 0.0)
      ! we need to wait for the ranks with ip_third = 0 - Barrier is implicit
      if (ip_third .eq. 0) then
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)
      endif
      ! We also need to broadcast halos
      call MPI_Bcast(ps, &
                     (ksizex_l + 2)*(ksizey_l + 2)*(ksizet_l + 2)*2, &
                     MPI_Real, 0, comm_grp_third, ierr)
      if (ip_third .eq. 0) then
        call MPI_WaitAll(12, reqs_pt, MPI_Statuses_Ignore, ierr)
      endif
      ! We also need to broadcast halos
      call MPI_Bcast(pt, &
                     (ksizex_l + 2)*(ksizey_l + 2)*(ksizet_l + 2)*2, &
                     MPI_Real, 0, comm_grp_third, ierr)
#else
      call gauss0(ps)
      psibarpsi1 = (0.0, 0.0)
      call gauss0(pt)
      psibarpsi2 = (0.0, 0.0)
#endif

      do idsource = 1, 2
        !
        !  source on domain wall at ithird=1
        !
        x = cmplx(0.0, 0.0)
        if (imass .ne. 5) then
          x(:, :, :, idsource) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))
        else
          x(:, :, :, idsource + 2) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))
        endif

        xi = cmplx(0.0, 0.0)
        if (ip_third .eq. 0) then
          if (imass .ne. 5) then
            xi(1, :, :, :, :) = x
          else
            !     xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
            xi(1, :, :, :, 1:2) = -x(:, :, :, 3:4)
            !     xi = 0.5(1+gamma_4)*eta on DW at ithird=1
          endif
        end if
#ifdef MPI
        call start_halo_update_5(4, xi, 11, reqs_xi)
        call complete_halo_update(reqs_xi)
#else
        call update_halo_5(4, xi)
#endif

        ! Phi= Mdagger*xi
        !
        call dslashd_shamir(Phi, xi, u, am, imass)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, Phi, 11, reqs_Phi)
        call complete_halo_update(reqs_Phi)
#else
        call update_halo_5(4, Phi)
#endif
        call congrad(Phi, res, itercg, am, imass)
        iter = iter + itercg

        if (imass .ne. 5) then
          if (ip_third .eq. np_third - 1) then
            !     pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
            psibarpsi1 = psibarpsi1 &
                         + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource)) &
                               *xi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))
          end if
        else
          if (ip_third .eq. 0) then
            !     pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
            psibarpsi1 = psibarpsi1 &
                         + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource + 2)) &
                               *xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource + 2))
          end if
        endif

        !  source on domain wall at ithird=kthird
        idsource2 = idsource + 2
        !
        x = cmplx(0.0, 0.0)

        if (imass .ne. 5) then
          x(:, :, :, idsource2) = cmplx(pt(:, :, :, 1), pt(:, :, :, 2))
        else
          x(:, :, :, idsource2 - 2) = cmplx(pt(:, :, :, 1), pt(:, :, :, 2))
        endif

        xi = cmplx(0.0, 0.0)
        if (ip_third .eq. np_third - 1) then
          if (imass .ne. 5) then
            !   xi = 0.5(1-gamma_4)*eta on DW at ithird=kthird
            xi(kthird_l, :, :, :, :) = x
          else
            !   xi = 0.5(1-gamma_4)*gamma_5*eta on DW at ithird=kthird
            xi(kthird_l, :, :, :, 3:4) = -x(:, :, :, 1:2)
          endif
        end if
#ifdef MPI
        call start_halo_update_5(4, xi, 11, reqs_xi)
        call complete_halo_update(reqs_xi)
#else
        call update_halo_5(4, xi)
#endif

        ! Phi= Mdagger*xi
        !
        call dslashd_shamir(Phi, xi, u, am, imass)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, Phi, 12, reqs_Phi)
        call complete_halo_update(reqs_Phi)
#else
        call update_halo_5(4, Phi)
#endif

        ! xi= (M)**-1 * Phi
        !
        call congrad(Phi, res, itercg, am, imass)
        iter = iter + itercg

        if (imass .ne. 5) then
          if (ip_third .eq. 0) then
            ! pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
            psibarpsi2 = psibarpsi2 &
                         + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2)) &
                               *xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2))
          end if
        else
          if (ip_third .eq. np_third - 1) then
            ! pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
            psibarpsi2 = psibarpsi2 &
                         + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource)) &
                               *xi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))
          end if
        endif

        !  end trace on Dirac indices....
      enddo

#ifdef MPI
      !  The sums psibarpsi1 and psibarpsi2 are initialised to 0 outside the loop
      !  So can be summed up per process, then collected here at the end
      call MPI_AllReduce(MPI_In_Place, psibarpsi1, 1, MPI_Double_Complex, &
                         MPI_Sum, comm, ierr)
      call MPI_AllReduce(MPI_In_Place, psibarpsi2, 1, MPI_Double_Complex, &
                         MPI_Sum, comm, ierr)
#endif

      if (imass .eq. 1) then
        psibarpsi1 = psibarpsi1/kvol
        psibarpsi2 = psibarpsi2/kvol
        pbp(inoise) = real(psibarpsi1 + psibarpsi2)
      elseif (imass .eq. 3) then
        psibarpsi1 = cmplx(0.0, -1.0)*psibarpsi1/kvol
        psibarpsi2 = cmplx(0.0, +1.0)*psibarpsi2/kvol
        pbp(inoise) = real(psibarpsi1 + psibarpsi2)
      elseif (imass .eq. 5) then
        psibarpsi1 = cmplx(0.0, -1.0)*psibarpsi1/kvol
        psibarpsi2 = cmplx(0.0, -1.0)*psibarpsi2/kvol
        pbp(inoise) = real(psibarpsi1 + psibarpsi2)
      endif

      !write(6,*) real(psibarpsi1),aimag(psibarpsi1), real(psibarpsi2),aimag(psibarpsi2)
      if (ip_global .eq. 0) then
        open (unit=100, file='fort.100', action='write', position='append')
        if (present(isweep_total)) then
          write (100, '(I5,4E17.9E3)') isweep_total, real(psibarpsi1), aimag(psibarpsi1), &
            real(psibarpsi2), aimag(psibarpsi2)
        else
          write (100, '(4E17.9E3)') real(psibarpsi1), aimag(psibarpsi1), &
            real(psibarpsi2), aimag(psibarpsi2)
        endif
        close (100)
      end if

      ! end loop on noise
    enddo

    susclsing = 0.0

    psibarpsi = real(sum(pbp))

    do inoise = 1, knoise
      susclsing = susclsing + real(sum(pbp(inoise)*pbp(inoise + 1:knoise)))
    enddo

    psibarpsi = psibarpsi/knoise
    susclsing = 2*kvol*susclsing/(knoise*(knoise - 1))

    if (ip_global .eq. 0) then
      open (unit=200, file='fort.200', action='write', position='append')
      if (present(isweep_total)) then
        write (200, '(I5,2E15.7E3)') isweep_total, psibarpsi, susclsing
      else
        write (200, '(2E15.7E3)') psibarpsi, susclsing
      endif
      close (200)
    end if

    aviter = float(iter)/(4*knoise)
    return

  end subroutine measure_shamir

  !*****************************************************************
  !   Calculate fermion expectation values via a noisy estimator
  !   -matrix inversion via conjugate gradient algorithm
  !       solves Mx=x1
  !     (Numerical Recipes section 2.10 pp.70-73)
  !*******************************************************************
  subroutine measure_wilson(psibarpsi, res, aviter, am, imass, isweep_total)
    use trial, only: u
    use vector, xi => x
    use comms5, only: start_halo_update_5
    use comms
    use gaussian
    use dirac
    use params
    implicit none
    real, intent(out) :: psibarpsi, aviter
    real, intent(in) :: res, am
    integer, intent(in) :: imass
    integer, intent(in), optional :: isweep_total
    complex(dp) :: x(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: psibarpsi1, psibarpsi2
    real(dp) :: cnum(0:1), cden(1)
    real :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
    real :: pt(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
    real(dp) :: pbp(knoise)
    integer :: idsource, idsource2, inoise
    integer :: iter, itercg
    real :: susclsing
#ifdef MPI
    integer, dimension(16) :: reqs_Phi, reqs_xi
    integer, dimension(12) :: reqs_ps, reqs_pt
    integer :: ierr
#endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Different from measure in master
    complex(dp) :: oslice(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: islice(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print *,"Wilson measure"

    iter = 0
    cnum(0) = 0.0
    cnum(1) = 1.0
    cden(1) = 0.0

    do inoise = 1, knoise
      !
      !     set up noise
#ifdef MPI
      if (ip_third .eq. 0) then
        call gauss0(ps, reqs_ps)
        call gauss0(pt, reqs_pt)
      endif
      psibarpsi1 = (0.0, 0.0)
      psibarpsi2 = (0.0, 0.0)
      ! we need to wait for the ranks with ip_third = 0 - Barrier is implicit
      if (ip_third .eq. 0) then
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)
      endif
      ! We also need to broadcast halos
      call MPI_Bcast(ps, &
                     (ksizex_l + 2)*(ksizey_l + 2)*(ksizet_l + 2)*2, &
                     MPI_Real, 0, comm_grp_third, ierr)
      if (ip_third .eq. 0) then
        call MPI_WaitAll(12, reqs_pt, MPI_Statuses_Ignore, ierr)
      endif
      ! We also need to broadcast halos
      call MPI_Bcast(pt, &
                     (ksizex_l + 2)*(ksizey_l + 2)*(ksizet_l + 2)*2, &
                     MPI_Real, 0, comm_grp_third, ierr)
#else
      call gauss0(ps)
      psibarpsi1 = (0.0, 0.0)
      call gauss0(pt)
      psibarpsi2 = (0.0, 0.0)
#endif

      do idsource = 1, 4
        !
        !  source on domain wall at ithird=1
        !
        x = cmplx(0.0, 0.0)
        x(:, :, :, idsource) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))

        xi = cmplx(0.0, 0.0)
        if (ip_third .eq. 0) then
          xi(1, :, :, :, :) = x
        end if
#ifdef MPI
        call start_halo_update_5(4, xi, 11, reqs_xi)
        call complete_halo_update(reqs_xi)
#else
        call update_halo_5(4, xi)
#endif
        !
        ! Phi= Mdagger*xi
        !
        call dslashd_wilson(Phi, xi, u, am, imass)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, Phi, 11, reqs_Phi)
        call complete_halo_update(reqs_Phi)
#else
        call update_halo_5(4, Phi)
#endif
        call congrad(Phi, res, itercg, am, imass)
        iter = iter + itercg
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This is not present in master
        ! for idirac=1 or 2
        islice=cmplx(0,0);
        islice(:,:,:,1:2)=xi(kthird_l,:,:,:,1:2) ! P+.phi
        call DWilson(oslice,islice,u,-am3)
        oslice=oslice-islice
        xi(kthird_l,:,:,:,:)=oslice
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        if (ip_third .eq. np_third - 1) then
          psibarpsi1 = psibarpsi1 &
                       + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource)) &
                             *xi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))
        end if
      end do 

      !  source on domain wall at ithird=kthird
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This is different from master 
      do idsource2 = 1, 4
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        x = cmplx(0.0, 0.0)

        x(:, :, :, idsource2) = cmplx(pt(:, :, :, 1), pt(:, :, :, 2))
        
        xi = cmplx(0.0, 0.0)
        if (ip_third .eq. np_third - 1) then
          xi(kthird_l, :, :, :, :) = x
        end if

#ifdef MPI
        call start_halo_update_5(4, xi, 11, reqs_xi)
        call complete_halo_update(reqs_xi)
#else
        call update_halo_5(4, xi)
#endif

        call dslashd_wilson(Phi, xi, u, am, imass)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, Phi, 12, reqs_Phi)
        call complete_halo_update(reqs_Phi)
#else
        call update_halo_5(4, Phi)
#endif

        ! xi= (M)**-1 * Phi
        !
        call congrad(Phi, res, itercg, am, imass)
        iter = iter + itercg

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This is not present in master
        ! for idirac=3 or 4
        islice=cmplx(0,0);
        islice(:,:,:,3:4)=xi(1,:,:,:,3:4) ! P-.phi
        call DWilson(oslice,islice,u,-am3)
        oslice=oslice-islice
        xi(1,:,:,:,:)=oslice
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (ip_third .eq. 0) then
          psibarpsi2 = psibarpsi2 &
                        + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2)) &
                              *xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2))
        end if

        !  end trace on Dirac indices....
      end do

#ifdef MPI
      !  The sums psibarpsi1 and psibarpsi2 are initialised to 0 outside the loop
      !  So can be summed up per process, then collected here at the end
      call MPI_AllReduce(MPI_In_Place, psibarpsi1, 1, MPI_Double_Complex, &
                         MPI_Sum, comm, ierr)
      call MPI_AllReduce(MPI_In_Place, psibarpsi2, 1, MPI_Double_Complex, &
                         MPI_Sum, comm, ierr)
#endif

      if (imass .eq. 1) then
        psibarpsi1 = psibarpsi1/kvol
        psibarpsi2 = psibarpsi2/kvol
        pbp(inoise) = real(psibarpsi1 + psibarpsi2)
      elseif (imass .eq. 3) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! The signs here are swapped compared to master 
        psibarpsi1 = cmplx(0.0, +1.0)*psibarpsi1/kvol
        psibarpsi2 = cmplx(0.0, -1.0)*psibarpsi2/kvol
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pbp(inoise) = real(psibarpsi1 + psibarpsi2)
      endif

      if (ip_global .eq. 0) then
        open (unit=100, file='fort.100', action='write', position='append')
        if (present(isweep_total)) then
          write (100, '(I5,4E17.9E3)') isweep_total, real(psibarpsi1), aimag(psibarpsi1), &
             & real(psibarpsi2), aimag(psibarpsi2)
        else
          write (100, '(4E17.9E3)') real(psibarpsi1), aimag(psibarpsi1), &
             & real(psibarpsi2), aimag(psibarpsi2)
        endif
        close (100)
      end if

      ! end loop on noise
    enddo

    susclsing = 0.0

    psibarpsi = real(sum(pbp))

    do inoise = 1, knoise
      susclsing = susclsing + real(sum(pbp(inoise)*pbp(inoise + 1:knoise)))
    enddo

    psibarpsi = psibarpsi/knoise
    susclsing = 2*kvol*susclsing/(knoise*(knoise - 1))

    if (ip_global .eq. 0) then
      open (unit=200, file='fort.200', action='write', position='append')
      if (present(isweep_total)) then
        write (200, '(I5,2E15.7E3)') isweep_total, psibarpsi, susclsing
      else
        write (200, '(2E15.7E3)') psibarpsi, susclsing
      endif
      close (200)
    end if

    aviter = float(iter)/(4*knoise)
    return

  end subroutine measure_wilson

  !******************************************************************
  !   Calculate meson correlators using point sources on domain walls
  !   -matrix inversion via conjugate gradient algorithm
  !       solves Mx=x1
  !     (Numerical Recipes section 2.10 pp.70-73)
  !*******************************************************************
  ! TODO: adjust calls to meson to supply output arrays
  ! DONE:: restored original call so that output written to disk from within
  ! subroutine SJH 3/11/21
  !subroutine meson(cpm, cmm, cferm1, cferm2, res, itercg, aviter, am, imass)
  ! output routines updated to be uniform with 'measure' subroutine
  !  SJH 7/12/21
  subroutine meson(res, itercg, aviter, am, imass, isweep_total)
    use random
    use vector, xi => x
    use dirac
    use trial
#ifdef MPI
    use comms4, only: start_halo_update_4
    use comms5, only: start_halo_update_5
#else
    use comms4, only: update_halo_4
    use comms5, only: update_halo_5

#endif

    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    real, intent(out) :: aviter
    integer, intent(in) :: imass
    integer, intent(in), optional :: isweep_total
    !! NOTICE : Full ksizet range.
    !real(dp), intent(out) :: cpm(0:ksizet - 1)
    !real(dp), intent(out) :: cmm(0:ksizet - 1)
    real(dp) :: cpm(0:ksizet - 1)
    real(dp) :: cmm(0:ksizet - 1)
    complex(dp) :: cpmn(0:ksizet - 1)
    complex(dp) :: cmmn(0:ksizet - 1)
    real(dp) :: tempcpmm_r(0:ksizet - 1)
    complex(dp) :: tempcpmm_c(0:ksizet - 1)
    !     complex x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
    !     complex xi,gamval
    !     complex prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
    complex(dp) :: x(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: x0(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: prop00(ksizex_l, ksizey_l, ksizet_l, 3:4, 1:2)
    complex(dp) :: prop0L(ksizex_l, ksizey_l, ksizet_l, 3:4, 3:4)
    complex(dp) :: prop00n(ksizex_l, ksizey_l, ksizet_l, 3:4, 1:2)
    complex(dp) :: prop0Ln(ksizex_l, ksizey_l, ksizet_l, 3:4, 3:4)
    !    complex :: cpmn(0:ksizet-1),cmmn(0:ksizet-1)
    !    real :: ps(ksizex_l, ksizey_l, ksizet_l, 2)
    real(dp) :: chim, chip
    integer :: ip_xxx, ip_yyy, ip_ttt, ixxx_l, iyyy_l, ittt_l
    integer :: ixxx, iyyy, ittt
    !     write(6,*) 'hi from meson'
    !
    integer, parameter :: nsource = 5
    integer, parameter :: nsmear = 0
    !    integer, parameter :: nsmear = 10
    real(dp), parameter :: c = 0.25d0
    integer :: iter, idsource, ksource, ismear
    integer :: it, itt, ittl
#ifdef MPI
    integer, dimension(16) :: mpireqs
    integer, dimension(12) :: mpireqs_4
    integer :: ierr
#endif
    !
    iter = 0
    itercg = 0
    !
    cpm = 0.0d0
    cmm = 0.0d0
    !
    cpmn = (0.0,0.0)
    cmmn = (0.0,0.0)
    !
    !      susceptibility
    chim = 0.0
    chip = 0.0
    !
    do ksource = 1, nsource
      !
      !   random location for +m source
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
#ifdef SITE_RANDOM
        ixxx = int(ksize*rano(yran, idum, 1, 1, 1)) + 1
        iyyy = int(ksize*rano(yran, idum, 1, 1, 1)) + 1
        ittt = int(ksizet*rano(yran, idum, 1, 1, 1)) + 1
#else
        ixxx = mod(2*isweep_total + ksource*(ksize/nsource), ksize) + 1
        iyyy = mod(3*isweep_total + ksource*(ksize/nsource), ksize) + 1
        ittt = mod(5*isweep_total + ksource*(ksizet/nsource), ksizet) + 1
#endif
#ifdef MPI
      endif
      call MPI_Bcast(ixxx, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(iyyy, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(ittt, 1, MPI_INTEGER, 0, comm, ierr)
#endif

      ip_xxx = int((ixxx - 1)/ksizex_l)
      ip_yyy = int((iyyy - 1)/ksizey_l)
      ip_ttt = int((ittt - 1)/ksizet_l)
      ixxx_l = mod(ixxx - 1, ksizex_l) + 1
      iyyy_l = mod(iyyy - 1, ksizey_l) + 1
      ittt_l = mod(ittt - 1, ksizet_l) + 1

      do idsource = 3, 4
        !  source on domain wall at ithird=1
        xi = (0.0d+0, 0.0d+0)
        x = (0.0d+0, 0.0d+0)
        !  point source at fixed site, spin...
        if (ip_xxx .eq. ip_x .and. ip_yyy .eq. ip_y .and. ip_ttt .eq. ip_t .and. ip_third .eq. 0) then
          x(ixxx_l, iyyy_l, ittt_l, idsource) = (1.0d+0, 0.0d+0)
        end if
        !
        ! now smear it.....
        !
#ifdef MPI
        call start_halo_update_4(4, x, 1, mpireqs_4)
        call MPI_Waitall(12, mpireqs_4, MPI_STATUSES_IGNORE, ierr)
#else
        call update_halo_4(4, x)
#endif

        do ismear = 1, nsmear
          call dslash2d(x0, x, u)
#ifdef MPI
          call start_halo_update_4(4, x0, 1, mpireqs_4)
          call MPI_Waitall(12, mpireqs_4, MPI_STATUSES_IGNORE, ierr)
#else
          call update_halo_4(4, x0)
#endif
          x = (1.0d+0 - c)*x + c*x0
        enddo! do ismear=1,nsmear
        !
        !
        !   xi = x  on DW at ithird=1
        !
        if (ip_third .eq. 0) then
          xi(1, :, :, :, :) = x
        endif
#ifdef MPI
        ! Overkill....
        call start_halo_update_5(4, xi, 1, mpireqs)
        call complete_halo_update(mpireqs)
#else
        ! Overkill....
        call update_halo_5(4, xi)
#endif

        !
        ! Phi= Mdagger*xi
        !
        call dslashd(Phi, xi, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, Phi, 1, mpireqs)
        call complete_halo_update(mpireqs)
#else
        call update_halo_5(4, Phi)
#endif
        !  preconditioning (no,really)
        call dslashd(xi, Phi, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, xi, 1, mpireqs)
        call complete_halo_update(mpireqs)
#else
        call update_halo_5(4, xi)
#endif
        !
        ! xi= (MdaggerM)**-1 * Phi
        !
        call congrad(Phi, res, itercg, am, imass)  ! solution is vector::x, here called xi
        iter = iter + itercg
        !
        if (ip_third .eq. 0) then
          prop00(:, :, :, idsource, 1:2) = xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
        else
          prop00(:, :, :, idsource, 1:2) = (0.0d0, 0.0d0)
        endif

        if (ip_third .eq. NP_THIRD - 1) then
          prop0L(:, :, :, idsource, 3:4) = xi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
        else
          prop0L(:, :, :, idsource, 3:4) = (0.0d0, 0.0d0)
        endif
        
        if (imass .ne. 1) then
          !  now evaluate with sign of mass reversed (not needed for hermitian mass term)
          !
          !   xi = x  on DW at ithird=1
          !
          xi = (0.0d+0, 0.0d+0)
          if (ip_third .eq. 0) then
            xi(1, :, :, :, :) = x
          endif
#ifdef MPI
          ! Overkill....
          call start_halo_update_5(4, xi, 1, mpireqs)
          call complete_halo_update(mpireqs)
#else
          ! Overkill....
          call update_halo_5(4, xi)
#endif

          !
          ! Phi= Mdagger*xi
          !
          call dslashd(Phi, xi, u, -am, imass)
#ifdef MPI
          call start_halo_update_5(4, Phi, 1, mpireqs)
          call complete_halo_update(mpireqs)
#else
          call update_halo_5(4, Phi)
#endif
          !  preconditioning (no,really)
          call dslashd(xi, Phi, u, -am, imass)
#ifdef MPI
          call start_halo_update_5(4, xi, 1, mpireqs)
          call complete_halo_update(mpireqs)
#else
          call update_halo_5(4, xi)
#endif
          !
          ! xi= (MdaggerM)**-1 * Phi
          !
          call congrad(Phi, res, itercg, -am, imass)  ! solution is vector::x, here called xi
          iter = iter + itercg
          !
          if (ip_third .eq. 0) then
            prop00n(:, :, :, idsource, 1:2) = xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
          else
            prop00n(:, :, :, idsource, 1:2) = (0.0d0, 0.0d0)
          endif

          if (ip_third .eq. NP_THIRD - 1) then
            prop0Ln(:, :, :, idsource, 3:4) = xi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
          else
            prop0Ln(:, :, :, idsource, 3:4) = (0.0d0, 0.0d0)
          endif

        else
          prop00n = prop00
          prop0Ln = prop0L
        endif
        
        !  end loop on source Dirac index....
      enddo ! do idsource=3,4


      ! Not actually necessary if result is used only by rank 0.
      call MPI_Bcast(prop00, size(prop00), MPI_DOUBLE_COMPLEX, &
                     0, comm_grp_third, ierr)

      call MPI_Barrier(comm_grp_third,ierr)

      call MPI_Bcast(prop00n, size(prop00n), MPI_DOUBLE_COMPLEX, &
                     0, comm_grp_third, ierr)

      call MPI_Barrier(comm_grp_third,ierr)

      call MPI_Bcast(prop0L, size(prop0L), MPI_DOUBLE_COMPLEX, &
                     NP_THIRD - 1, comm_grp_third, ierr)

      call MPI_Barrier(comm_grp_third,ierr)

      call MPI_Bcast(prop0Ln, size(prop0Ln), MPI_DOUBLE_COMPLEX, &
                     NP_THIRD - 1, comm_grp_third, ierr)

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !
      !  Now tie up the ends....
      !
      !  First C+-
      !
      !  now evaluate the trace (exploiting projection)
      tempcpmm_r = 0.d0
      tempcpmm_c = (0.d0, 0.d0)
      do it = 0, ksizet - 1
        itt = mod((ittt + it - 1), ksizet) + 1
        ittl = itt - ip_t*ksizet_l
        if (ittl .ge. 1 .and. ittl .le. ksizet_l) then
          tempcpmm_r(it) = sum(abs(prop00(:, :, ittl, 3:4, 1:2))**2)
          tempcpmm_c(it) = sum(prop00(:, :, ittl, 3:4, 1:2)*conjg(prop00n(:, :, ittl, 3:4, 1:2)))
        endif
      enddo
#ifdef MPI
      timered1: block
        integer :: comm_grp_third_dual
        call MPI_Comm_split(MPI_COMM_WORLD, ip_third, &
                            ip_x + ip_y*np_x + ip_t*np_x*np_y, comm_grp_third_dual, ierr)

        call MPI_AllReduce(MPI_In_Place, tempcpmm_r, ksizet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_grp_third_dual, ierr)
        call MPI_AllReduce(MPI_In_Place, tempcpmm_c, ksizet, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_grp_third_dual, ierr)
      end block timered1
      if (ip_global .eq. 0) then
#endif
        cpm = cpm + tempcpmm_r
        cpmn = cpmn + tempcpmm_c
#ifdef MPI
      endif! if(ip_global.eq.0) then
#endif
      !
      !  next C--
      !  now evaluate the trace exploiting projection
      tempcpmm_r = 0.d0
      tempcpmm_c = (0.d0, 0.d0)
      do it = 0, ksizet - 1
        itt = mod((ittt + it - 1), ksizet) + 1
        ittl = itt - ip_t*ksizet_l
        if (ittl .ge. 1 .and. ittl .le. ksizet_l) then
          tempcpmm_r(it) = sum(abs(prop0L(:, :, ittl, 3:4, 3:4))**2)
          tempcpmm_c(it) = sum(prop0L(:, :, ittl, 3:4, 3:4)*conjg(prop0Ln(:, :, ittl, 3:4, 3:4)))
        endif
      enddo
#ifdef MPI
      timered2: block
        integer :: comm_grp_third_dual
        call MPI_Comm_split(MPI_COMM_WORLD, ip_third, &
                            ip_x + ip_y*np_x + ip_t*np_x*np_y, comm_grp_third_dual, ierr)

        call MPI_AllReduce(MPI_In_Place, tempcpmm_r, ksizet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_grp_third_dual, ierr)
        call MPI_AllReduce(MPI_In_Place, tempcpmm_c, ksizet, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_grp_third_dual, ierr)
      end block timered2

      if (ip_global .eq. 0) then
#endif
        cmm = cmm + tempcpmm_r
        cmmn = cmmn + tempcpmm_c
#ifdef MPI
      endif
#endif

      !  finish loop over sources
    enddo! do ksource=1,nsource
    !
    do it = 0, ksizet - 1
      cpm(it) = cpm(it)/nsource
      cmm(it) = cmm(it)/nsource
      !  Cf. (54) of 1507.07717
      cpmn(it) = cpmn(it)/nsource
      cmmn(it) = cmmn(it)/nsource
      !  Cf. (60) of 1507.07717
      chim = chim + 2*(cpm(it) + cmm(it))
    enddo

#ifdef MPI
    if (ip_global .eq. 0) then
#endif
        open (unit=302, file='fort.302', action='write', position='append')
        if (present(isweep_total)) then
          do it = 0, ksizet - 1
             write (302, *) isweep_total, it, cpm(it), cmm(it)
          enddo
        else
          do it = 0, ksizet-1
             write (302, *) it, cpm(it), cmm(it)
          enddo
        endif
        close (302)
        open (unit=320, file='fort.320', action='write', position='append')
        if (present(isweep_total)) then
          do it = 0, ksizet - 1
             write (320, *) isweep_total, it, real(cpmn(it)), aimag(cpmn(it))
          enddo
        else
          do it = 0, ksizet-1
             write (320, *) it, real(cpmn(it)), aimag(cpmn(it))
          enddo
        endif
        close (320)
        open (unit=321, file='fort.321', action='write', position='append')
        if (present(isweep_total)) then
          do it = 0, ksizet - 1
             write (321, *) isweep_total, it, real(cmmn(it)), aimag(cmmn(it))
          enddo
        else
          do it = 0, ksizet-1
             write (321, *) it, real(cmmn(it)), aimag(cmmn(it))
          enddo
        endif
        close (321)
        open (unit=400, file='fort.400', action='write', position='append')
        if (present(isweep_total)) then
             write (400, *) isweep_total, chim
        else
             write (400, *) chim
        endif
        close (400)
      !
      !     write(6,*) chim
#ifdef MPI
    endif ! if(ip_global.eq.0) then
#endif
    !     if(imass.ne.1)then
    !     do it=0,ksizet-1
    !     write(402,*) it, real(cpmn(it)), real(cmmn(it))
    !     write(403,*) it, aimag(cpmn(it)), aimag(cmmn(it))
    !     enddo
    !     write(401,*) chip
    !     endif
    !
    !     if(imass.eq.1)then
    aviter = float(iter)/(2*nsource)
    !     else
    !     aviter=float(iter)/(4*nsource)
    !     endif
    !
    return
  end subroutine meson
!******************************************************************

  !******************************************************************
  !   Calculate fermion correlators using wall sources on domain walls
  !   -matrix inversion via conjugate gradient algorithm
  !       solves Mx=x1
  !     (Numerical Recipes section 2.10 pp.70-73)
  !*******************************************************************
  ! TODO: adjust calls to fermion to supply output arrays
  ! DONE:: restored original call so that output written to disk from within
  ! subroutine SJH 3/11/21
  !subroutine fermion(cpm, cmm, cferm1, cferm2, res, itercg, aviter, am, imass)
  ! output routines updated to be uniform with 'measure' subroutine
  !  SJH 7/12/21
  subroutine fermion(res, itercg, aviter, am, imass, isweep_total)
    use random
    use vector, xi => x
    use dirac
    use trial
#ifdef MPI
    use comms4, only: start_halo_update_4
    use comms5, only: start_halo_update_5
#else
    use comms4, only: update_halo_4
    use comms5, only: update_halo_5

#endif

    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    real, intent(out) :: aviter
    integer, intent(in) :: imass
    integer, intent(in), optional :: isweep_total
    !! NOTICE : Full ksizet range.
    complex(dp) :: tempcpmm_c(0:ksizet - 1)
    !complex(dp), intent(out) :: cferm1(0:ksizet - 1)
    !complex(dp), intent(out) :: cferm2(0:ksizet - 1)
    complex(dp) :: cferm1(0:ksizet - 1)
    complex(dp) :: cferm2(0:ksizet - 1)
    !     complex xi,gamval
    !     complex prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
    complex(dp) :: x(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: prop00(ksizex_l, ksizey_l, ksizet_l, 3:4, 1:2)
    complex(dp) :: prop0L(ksizex_l, ksizey_l, ksizet_l, 3:4, 3:4)
    !    complex :: prop00n(ksizex_l, ksizey_l, ksizet_l, 3:4, 1:2)
    !    complex :: prop0Ln(ksizex_l, ksizey_l, ksizet_l, 3:4, 3:4)
    !    complex :: cpmn(0:ksizet-1),cmmn(0:ksizet-1)
    !    real :: ps(ksizex_l, ksizey_l, ksizet_l, 2)
    real(dp) :: chim, chip
    integer :: ip_ttt, ittt_l
    integer :: ittt
    !     write(6,*) 'hi from fermion'
    !
    integer, parameter :: nsource = 5
    real(dp), parameter :: c = 0.25d0
    integer :: iter, idsource, ksource, isign
    integer :: it, itt, ittl, idd
#ifdef MPI
    integer, dimension(16) :: mpireqs
    integer, dimension(12) :: mpireqs_4
    integer :: ierr
#endif
    !
    iter = 0
    itercg = 0
    !
    cferm1 = (0.0d+0, 0.0d+0)
    cferm2 = (0.0d+0, 0.0d+0)
    !
    !      susceptibility
    chim = 0.0
    chip = 0.0
    !
    do ksource = 1, nsource
      !
      !   random location for +m source
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
#ifdef SITE_RANDOM
        ittt = int(ksizet*rano(yran, idum, 1, 1, 1)) + 1
#else
        ittt = mod(5*isweep_total + ksource*(ksizet/nsource), ksizet) + 1
#endif
#ifdef MPI
      endif
      call MPI_Bcast(ittt, 1, MPI_INTEGER, 0, comm, ierr)
#endif

      ip_ttt = int((ittt - 1)/ksizet_l)
      ittt_l = mod(ittt - 1, ksizet_l) + 1

      do idsource = 3, 4
        !  source on domain wall at ithird=1
        xi = (0.0d+0, 0.0d+0)
        x = (0.0d+0, 0.0d+0)
        !  wall source
        if (ip_ttt .eq. ip_t .and. ip_third .eq. 0) then
          x(:, :, ittt_l, idsource) = cmplx(1.0,0.0) / ksize2
        end if

#ifdef MPI
        call start_halo_update_4(4, x, 1, mpireqs_4)
        call MPI_Waitall(12, mpireqs_4, MPI_STATUSES_IGNORE, ierr)
#else
        call update_halo_4(4, x)
#endif
        !
        !   xi = x  on DW at ithird=1
        !
        if (ip_third .eq. 0) then
          xi(1, :, :, :, :) = x
        endif
#ifdef MPI
        ! Overkill....
        call start_halo_update_5(4, xi, 1, mpireqs)
        call complete_halo_update(mpireqs)
#else
        ! Overkill....
        call update_halo_5(4, xi)
#endif

        !
        ! Phi= Mdagger*xi
        !
        call dslashd(Phi, xi, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, Phi, 1, mpireqs)
        call complete_halo_update(mpireqs)
#else
        call update_halo_5(4, Phi)
#endif
        !  preconditioning (no,really)
        call dslashd(xi, Phi, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, xi, 1, mpireqs)
        call complete_halo_update(mpireqs)
#else
        call update_halo_5(4, xi)
#endif
        !
        ! xi= (MdaggerM)**-1 * Phi
        !
        call congrad(Phi, res, itercg, am, imass)  ! solution is vector::x, here called xi
        iter = iter + itercg
        !
        if (ip_third .eq. 0) then
          prop00(:, :, :, idsource, 1:2) = xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
        else
          prop00(:, :, :, idsource, 1:2) = (0.0d0, 0.0d0)
        endif

        if (ip_third .eq. NP_THIRD - 1) then
          prop0L(:, :, :, idsource, 3:4) = xi(kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
        else
          prop0L(:, :, :, idsource, 3:4) = (0.0d0, 0.0d0)
        endif

        !  end loop on source Dirac index....
      enddo ! do idsource=3,4


      ! Not actually necessary if result is used only by rank 0.
      call MPI_Bcast(prop00, size(prop00), MPI_DOUBLE_COMPLEX, &
                     0, comm_grp_third, ierr)

      call MPI_Barrier(comm_grp_third,ierr)

      call MPI_Bcast(prop0L, size(prop0L), MPI_DOUBLE_COMPLEX, &
                     NP_THIRD - 1, comm_grp_third, ierr)

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !
      !  Now tie up the ends....
      !
      !    now the fermion propagator
      !  = tr{ P_-*Psi(0,1)Psibar(x,Ls) + gamma_0*P_-*Psi(0,1)Psibar(x,1) }
      tempcpmm_c = (0.0d+0, 0.0d+0)
      do idd = 3, 4
        do it = 0, ksizet - 1
          itt = mod((ittt + it - 1), ksizet) + 1
          ! correct for apbc
          if (itt .ge. ittt) then
            isign = 1
          else
            isign = ibound
          endif
          ittl = itt - ip_t*ksizet_l
          ! Working only on the local T range
          if (ittl .ge. 1 .and. ittl .le. ksizet_l) then
            tempcpmm_c(it) = tempcpmm_c(it) + &
              & isign*akappa*sum(prop0L(:, :, ittl, idd, idd))
          endif
        enddo
      enddo! do idd=3,4
#ifdef MPI
      timered3: block
        integer :: comm_grp_third_dual
        call MPI_Comm_split(MPI_COMM_WORLD, ip_third, &
                            ip_x + ip_y*np_x + ip_t*np_x*np_y, comm_grp_third_dual, ierr)

        call MPI_AllReduce(MPI_In_Place, tempcpmm_c, ksizet, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_grp_third_dual, ierr)
      end block timered3

      if (ip_global .eq. 0) then
#endif
        !!! if(ip_global.eq.0) then
        cferm1 = cferm1 + tempcpmm_c
#ifdef MPI
      endif
#endif
      tempcpmm_c = (0.0d0, 0.0d0)
      do idd = 3, 4
        do it = 0, ksizet - 1
          itt = mod((ittt + it - 1), ksizet) + 1
          ! correct for apbc
          if (itt .ge. ittt) then
            isign = 1
          else
            isign = ibound
          endif
          ittl = itt - ip_t*ksizet_l
          ! Working only on the local T range
          if (ittl .ge. 1 .and. ittl .le. ksizet_l) then
            tempcpmm_c(it) = tempcpmm_c(it) + &
             & isign*gamval(3, idd)*sum(prop00(:, :, ittl, idd, gamin(3, idd)))
          endif
        enddo
      enddo! do idd=3,4
#ifdef MPI
      timered4: block
        integer :: comm_grp_third_dual
        call MPI_Comm_split(MPI_COMM_WORLD, ip_third, &
                            ip_x + ip_y*np_x + ip_t*np_x*np_y, comm_grp_third_dual, ierr)

        call MPI_AllReduce(MPI_In_Place, tempcpmm_c, ksizet, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_grp_third_dual, ierr)
      end block timered4

      if (ip_global .eq. 0) then
#endif
        !!! if(ip_global.eq.0) then
        cferm2 = cferm2 + tempcpmm_c
#ifdef MPI
      endif
#endif

      !  finish loop over sources
    enddo! do ksource=1,nsource
    !
    do it = 0, ksizet - 1
      cferm1(it) = cferm1(it)/nsource
      cferm2(it) = cferm2(it)/nsource
    enddo

#ifdef MPI
    if (ip_global .eq. 0) then
#endif
        open (unit=500, file='fort.500', action='write', position='append')
        if (present(isweep_total)) then
          do it = 0, ksizet - 1
             write (500, *) isweep_total, it, real(cferm1(it)), aimag(cferm1(it))
          enddo
        else
          do it = 0, ksizet - 1
             write (500, *) it, real(cferm1(it)), aimag(cferm1(it))
          enddo
        endif
        close (500)
        open (unit=501, file='fort.501', action='write', position='append')
        if (present(isweep_total)) then
          do it = 0, ksizet - 1
             write (501, *) isweep_total, it, real(cferm2(it)), aimag(cferm2(it))
          enddo
        else
          do it = 0, ksizet - 1
             write (501, *) it, real(cferm2(it)), aimag(cferm2(it))
          enddo
        endif
        close (501)
      !
      !     write(6,*) chim
#ifdef MPI
    endif ! if(ip_global.eq.0) then
#endif

    aviter = float(iter)/(2*nsource)

    return

  end subroutine fermion
!******************************************************************

end module measure_module
