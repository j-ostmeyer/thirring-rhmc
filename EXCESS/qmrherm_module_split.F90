module qmrherm_module_split
  use params
  implicit none
  complex(dp) :: vtild(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: q(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: pm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: qm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: p(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x3(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: x1(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

  complex(dp), save :: Phi0(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, ndiag)
  logical :: printall

contains
  !******************************************************************
  !    multisolver matrix inversion via Lanczos technique
  !  eg. Golub & van Loan "Matrix Computations" 9.3.1
  ! http://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf
  !       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
  !   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
  !   can be interchanged with congrad for p=-1
  !   iflag=1: in addition updates Phi0 register needed for PV force term
  !   iflag=2: evaluates DWF force term
  !   iflag=3: evaluates PV force term
  !*****************************************************************m
  subroutine qmrherm_split(Phi, X, res, itercg, am, imass, anum, aden, ndiagq, iflag, isweep, &
      & iter)
    use params
    use trial, only: u
    use gforce
    use comms, only: complete_halo_update
    use comms5, only: start_halo_update_5
    use comms6, only: start_halo_update_6
    use partitioning
    use comms_partitioning
    use dirac_split
    use dirac
    use derivs_module
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(out) :: X(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: imass, ndiagq, iflag, isweep, iter
    real(dp), intent(in) :: anum(0:ndiagq), aden(ndiagq)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    !
    real(dp) :: alphatild
    real(dp) :: coeff
    !
    real(dp) :: alpha(ndiagq)
    real(dp) :: amu(ndiagq), d(ndiagq), dm1(ndiagq)
    real(dp) :: rho(ndiagq), rhom1(ndiagq)
    real(dp) :: betaq, betaq0, phimod
    real :: resid, rhomax, arelax
    integer :: niter, idiag
    logical :: go_on
#ifdef MPI
    integer, dimension(12) :: reqs_X2, reqs_Phi0, reqs_R, reqs_x
    integer :: ierr
    real(dp) :: dp_reduction ! DEBUG
    ! For dirac split
    integer :: Rsreqs(54), Rrreqs(54)
    integer :: vtildsreqs(54), vtildrreqs(54)
#endif
    integer :: ichunk(3), mu
    integer :: wpc
    ! initialize communications
    call init_partitioning
    ! MPI datatypes for send and receive
    call init_dirac_hb_types
    ! MPI_requests
    call get_dirac_recvreqs(Rrreqs, R)
    call get_dirac_sendreqs(Rsreqs, R)
    call get_dirac_recvreqs(vtildrreqs, vtild)
    call get_dirac_sendreqs(vtildsreqs, vtild)
    ! list of work
    call get_dslash_work_ordering(dslash_work_ordering, .false.)
    call get_dslash_work_ordering(dslashd_work_ordering, .true.)

    resid = sqrt(kthird*ksize*ksize*ksizet*4*res*res)

    itercg = 0
    !
    !   initialise r=Phi
    !
    R = Phi
    qm1 = cmplx(0.0, 0.0)
    x = anum(0)*Phi
    betaq = sum(abs(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)

    !call MPI_AllReduce(MPI_In_Place, betaq, 1, MPI_Double_Precision, MPI_Sum, comm,ierr) ! DEBUG
    call MPI_AllReduce(betaq, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr) ! DEBUG
    betaq = dp_reduction

    betaq = sqrt(betaq)
    phimod = betaq
    !
    !do niter=1,20
    niter = 0
    go_on = .true.
    call MPI_StartAll(54, Rrreqs, ierr)
    call MPI_StartAll(54, Rsreqs, ierr)

    do while (niter .lt. max_qmr_iters .and. go_on)
      niter = niter + 1
      itercg = itercg + 1
      !
      !  Lanczos steps
      !
      ! on the local lattice
      q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
        & R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)/betaq
      ! q = R / betaq on the halo && call dslash(vtild,q,u,am,imass)
      ! starting recv requests for vtild
      ! send requests are started in hbetaqdiv_dslash_split
      call MPI_StartAll(54, vtildrreqs, ierr)
      dslash_swd = .false.
      do wpc = 1, 27*7
        ichunk = dslash_work_ordering(1:3, wpc)
        mu = dslash_work_ordering(4, wpc)
        call hbetaqdiv_dslash_split(q, betaq, vtild, R, u, am, imass, ichunk, mu,&
             & border_partitions_cube, dslash_swd, Rrreqs, vtildsreqs)
      enddo
      call MPI_Barrier(comm, ierr)

      alphatild = sum(abs(vtild(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)**2))
      !call MPI_AllReduce(MPI_In_Place, alphatild, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
      call MPI_AllReduce(alphatild, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
      alphatild = dp_reduction

      ! starting recv requests for R
      ! send requests are started in dslashd_Rcomp_split
      call MPI_StartAll(54, Rrreqs, ierr)
      dslashd_swd = .false.
      do wpc = 1, 27*7
        ichunk = dslashd_work_ordering(1:3, wpc)
        mu = dslashd_work_ordering(4, wpc)
        call dslashd_Rcomp_split(R, x3, alphatild, q, betaq, qm1, vtild, u, am, imass, ichunk, mu,&
          & border_partitions_cube, dslashd_swd, vtildrreqs, Rsreqs)
      enddo

      call MPI_Barrier(comm, ierr)

      qm1 = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
      !
      betaq0 = betaq

      betaq = sum(abs(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
      !call MPI_AllReduce(MPI_In_Place, betaq, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
      call MPI_AllReduce(betaq, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
      betaq = dp_reduction
      betaq = sqrt(betaq)
      !print*,'betaq',betaq,ip_global
      !
      alpha = alphatild + aden
      !print*,'alpha',alpha,ip_global
      call MPI_Barrier(comm, ierr)

      !
      if (niter .eq. 1) then
        d = alpha
        rho = betaq0/alpha
        rhom1 = rho
        do idiag = 1, ndiagq
          p(:, :, :, :, :, idiag) = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
          x1(:, :, :, :, :, idiag) = rho(idiag)*p(:, :, :, :, :, idiag)
        enddo
        pm1 = p
      else
        amu = betaq0/d
        dm1 = d
        d = alpha - betaq0*amu
        rho = -amu*dm1*rhom1/d
        do idiag = 1, ndiagq
          p(:, :, :, :, :, idiag) = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
              & - amu(idiag)*pm1(:, :, :, :, :, idiag)
        enddo
        pm1 = p
        !     Convergence criterion (a bit ad hoc for now...)
        rhomax = real(maxval(abs(phimod*rho)))
        rhom1 = rho
        do idiag = 1, ndiagq
          x1(:, :, :, :, :, idiag) = &
            & x1(:, :, :, :, :, idiag) &
            & + rho(idiag)*p(:, :, :, :, :, idiag)
        enddo

        !     check to see whether the residual is acceptable for all ndiagq....
        !     criterion is a bit ad hoc -- relaxing by a factor arelax improves code
        !     stability and leads to quicker convergence
        arelax = 2.0
        if (rhomax .lt. arelax*resid) then
          !     if(rhomax.lt.resid) then
          !     call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiagq)
          !     convergence based on || residual || not working well in single precision...
          !     if(resmax.lt.resid) goto 8
          go_on = .false.
        endif
      endif
      !
    enddo! do while(niter.lt.max_qmr_iters .and. go_on )

    if (niter .gt. max_qmr_iters) then
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
        write (7, *) 'QMRniterc!, niter, isweep,iter,iflag,imass,anum,ndiagq = ', &
            &   niter, isweep, iter, iflag, imass, anum(0), ndiagq
#ifdef MPI
      end if
#endif
    endif
    !
    !8   continue
    if (iflag .lt. 2) then
      !     Now evaluate solution x=(MdaggerM)^p * Phi
      do idiag = 1, ndiagq
        x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
          & x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
          & + anum(idiag)*x1(:, :, :, :, :, idiag)
      enddo
#ifdef MPI
      ! x is a saved module variable, so must be updated to avoid polluting the parent function
      ! could this in principle be moved outside the function so we don't do it unnecessarily?
      ! but in that case we wouldn't be able to hide the communications
      call start_halo_update_5(4, x, 3, reqs_x)
#else
      call update_halo_5(4, x)
#endif
      !
      !  update phi0 block if required...
      if (iflag .eq. 1) then
        Phi0(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, 1:ndiagq) = &
          & X1(:, :, :, :, :, 1:ndiagq)
#ifdef MPI
        ! No way to hide communications here unfortunately
        ! In principle this could be better interleaved with the x update
        ! but that would add extra branching, and this section is messy enough already
        call start_halo_update_6(4, ndiagq, Phi0, 4, reqs_Phi0)
        call complete_halo_update(reqs_Phi0)
#else
        call update_halo_6(4, ndiagq, Phi0)
#endif
      endif! if(iflag.eq.1) then
#ifdef MPI
      call complete_halo_update(reqs_x)
#endif
      !
    else! if(iflag.lt.2)then
      !
      do idiag = 1, ndiagq
        !
        !  X2 = M*X1
        R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = X1(:, :, :, :, :, idiag)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, R, 5, reqs_R)
        call complete_halo_update(reqs_R)
#else
        call update_halo_5(4, R)
#endif

        ! Communication of X2 generated here can be hidden if iflag isn't 2, while R is updated
        call dslash(X2, R, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, X2, 6, reqs_X2)
#else
        call update_halo_5(4, X2)
#endif
        !
        if (iflag .eq. 2) then
          coeff = anum(idiag)
#ifdef MPI
          call complete_halo_update(reqs_X2)
#endif
          call derivs(R, X2, coeff, 0)
        else! if(iflag.eq.2)then
          coeff = -anum(idiag)
          R = Phi0(:, :, :, :, :, idiag)
#ifdef MPI
          call complete_halo_update(reqs_X2)
#endif
          call derivs(R, X2, coeff, 0)
          !
          ! Communication of X2 generated here can be hidden while R is updated
          call dslash(X2, R, u, am, imass)
#ifdef MPI
          call start_halo_update_5(4, X2, 7, reqs_X2)
#else
          call update_halo_5(4, X2)
#endif
          !
          R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(:, :, :, :, :, idiag)
#ifdef MPI
          call start_halo_update_5(4, R, 8, reqs_R)
          call complete_halo_update(reqs_X2)
          call complete_halo_update(reqs_R)
#else
          call update_halo_5(4, R)
#endif
          call derivs(X2, R, coeff, 1)
        endif! if(iflag.eq.2)then
      enddo! do idiag=1, ndiagq
    endif !if(iflag.lt.2)then , else

    if (ip_global .eq. 0 .and. printall) then
      print *, "Qmrherm split iterations,res:", itercg, res
    endif
    return
  end subroutine qmrherm_split
  !**********************************************************************
  !  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
  !  iflag = 1 : evaluates Rdagger*(M)'*X2
  !**********************************************************************

  ! this subroutine merges the q = R / betaq step and dslash application
  ! WARNING: the q = R/ betaq is performed only on the halo,
  !          it is assumed that q is already computed on the local lattice.
  !          This complication arises from the fact that the values of q needed
  !          for the dslash_calculation can be in a different partition from
  !          the one specified by ichunk - shifted by +-1 site in one of the 3
  !          directions.

  subroutine hbetaqdiv_dslash_split(tq, betaq, Phi, R, u, am, imass, ichunk, mu, tbpc, tdsswd, tdhrr, tdbsr)
    use params
    use partitioning
    use mpi
    use dirac_split
    implicit none
    complex(dp), intent(inout) :: tq(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: betaq
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: ichunk(3) ! portion of array to operate on
    integer, intent(in) :: mu ! -3 <= mu <= 3
    ! Temp Border Partition Cube
    type(localpart), intent(in) :: tbpc(-1:1, -1:1, -1:1)
    ! Temp DSlash Split Work Done
    logical, intent(inout) :: tdsswd(-3:3, -1:1, -1:1, -1:1)
    ! Temp Dirac Halo Recv Requests
    integer, intent(inout) :: tdhrr(54)
    integer, intent(inout) :: tdbsr(54)

    integer :: chunk(2, 3)
    ! CHUNK Shifted
    integer :: chunk_s(2, 3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: inn
    integer :: xd, xu, yd, yu, td, tu ! portion of array to operate on (Phi)
    integer :: ierr

    tpart = tbpc(ichunk(1), ichunk(2), ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not. any(tdsswd(:, ichunk(1), ichunk(2), ichunk(3)))

    if (halo_to_wait_for .ne. 0) then
      call MPI_Wait(tdhrr(halo_to_wait_for), MPI_STATUS_IGNORE, ierr)

      ! performing q = R/betaq first
      chunk_s = chunk
      chunk_s(:, abs(mu)) = chunk_s(:, abs(mu)) + sign(1, mu)
      xd = chunk_s(1, 1)
      xu = chunk_s(2, 1)
      yd = chunk_s(1, 2)
      yu = chunk_s(2, 2)
      td = chunk_s(1, 3)
      tu = chunk_s(2, 3)

      tq(:, xd:xu, yd:yu, td:tu, :) = R(:, xd:xu, yd:yu, td:tu, :)/betaq
    endif

    if (mu .eq. 0) then
      call dslash_split_local(Phi, tq, am, imass, chunk, init)
    else
      if (mu .gt. 0) then
        call dslash_split_nonlocal(Phi, tq, u, chunk, mu, 1, init)
      else if (mu .lt. 0) then
        call dslash_split_nonlocal(Phi, tq, u, chunk, -mu, -1, init)
      endif
    endif

    ! flagging work done
    tdsswd(mu, ichunk(1), ichunk(2), ichunk(3)) = .true.
    if (all(tdsswd(:, ichunk(1), ichunk(2), ichunk(3)))) then
      tpart = tbpc(ichunk(1), ichunk(2), ichunk(3))
      do inn = 1, tpart%nn
        ! clearing send requests
        call MPI_Wait(tdbsr(tpart%ahpss(inn)), MPI_STATUS_IGNORE, ierr)
        ! restarting send request
        call MPI_Start(tdbsr(tpart%ahpss(inn)), ierr)
      enddo
    endif

  end subroutine

  ! this subroutine merges the dslashd application and the
  ! R = x3 - alphatild * q - betaq * qm1
  ! steps
  subroutine dslashd_Rcomp_split(R, x3, alphatild, tq, betaq, tqm1, vtild, u, am, imass,&
                          & ichunk, mu, tbpc, tdsswd, tdhrr, tdbsr)
    use params
    use partitioning
    use mpi
    use dirac_split
    use comms ! DEBUG
    implicit none
    complex(dp), intent(out) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(out) :: x3(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: alphatild
    complex(dp), intent(in) :: tq(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp), intent(in) :: betaq
    complex(dp), intent(in) :: tqm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp), intent(in) :: vtild(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: ichunk(3) ! portion of array to operate on
    integer, intent(in) :: mu ! -3 <= mu <= 3
    ! Temp Border Partition Cube
    type(localpart), intent(in) :: tbpc(-1:1, -1:1, -1:1)
    ! Temp DSlash Split Work Done
    logical, intent(inout) :: tdsswd(-3:3, -1:1, -1:1, -1:1)
    ! Temp Dirac Halo Recv Requests
    integer, intent(inout) :: tdhrr(54)
    ! Temp Dirac Border Send Requests
    integer, intent(inout) :: tdbsr(54)

    integer :: chunk(2, 3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: inn
    integer :: xd, xu, yd, yu, td, tu ! portion of array to operate on
    integer :: ierr

    tpart = tbpc(ichunk(1), ichunk(2), ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not. any(tdsswd(:, ichunk(1), ichunk(2), ichunk(3)))

    if (halo_to_wait_for .ne. 0) then
      call MPI_Wait(tdhrr(halo_to_wait_for), MPI_STATUS_IGNORE, ierr)
    endif
    if (mu .eq. 0) then
      call dslashd_split_local(x3, vtild, am, imass, chunk, init)
    else
      if (mu .gt. 0) then
        call dslashd_split_nonlocal(x3, vtild, u, chunk, mu, 1, init)
      else if (mu .lt. 0) then
        call dslashd_split_nonlocal(x3, vtild, u, chunk, -mu, -1, init)
      endif
    endif

    ! R = x3 - alphatild * q - betaq * qm1
    ! on the required partition
    xd = chunk(1, 1)
    xu = chunk(2, 1)
    yd = chunk(1, 2)
    yu = chunk(2, 2)
    td = chunk(1, 3)
    tu = chunk(2, 3)

    R(:, xd:xu, yd:yu, td:tu, :) = x3(:, xd:xu, yd:yu, td:tu, :) - &
      & alphatild*tq(:, xd:xu, yd:yu, td:tu, :) - &
      & betaq*tqm1(:, xd:xu, yd:yu, td:tu, :)

    ! flagging work done
    tdsswd(mu, ichunk(1), ichunk(2), ichunk(3)) = .true.
    ! checking whether to send the partition already or not
    if (all(tdsswd(:, ichunk(1), ichunk(2), ichunk(3)))) then
      tpart = tbpc(ichunk(1), ichunk(2), ichunk(3))
      do inn = 1, tpart%nn
        ! clearing send requests
        call MPI_Wait(tdbsr(tpart%ahpss(inn)), MPI_STATUS_IGNORE, ierr)
        ! restarting send request
        call MPI_Start(tdbsr(tpart%ahpss(inn)), ierr)
      enddo
    endif

  end subroutine

end module qmrherm_module_split
