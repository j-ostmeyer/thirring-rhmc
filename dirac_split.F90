module dirac_split
! Versions of the dirac operator subroutines which operate only on a subset 
! of the lattice. 
  use params
  use comms
  implicit none
  save

  logical :: dslash_swd(-3:3,-1:1,-1:1,-1:1)  ! DSLASH Split Work Done
  logical :: dslashd_swd(-3:3,-1:1,-1:1,-1:1) ! DSLASHD Split Work Done

  ! list of partitions, indexed by ipx,ipy and ipt partition coordinates, 
  ! and mu.
  ! with  -1 <= ip[xyt] <= 1, -3 <= mu <= 3
  ! in principle it can also depend on the subroutine they are used in 
  ! (qmrherm or congrad)
  integer :: dslash_work_ordering(4,27*7)
  integer :: dslashd_work_ordering(4,27*7)

contains 

! DSLASH 

  subroutine dslash_split(Phi,R,u,am,imass,ichunk,mu,tbpc,tdsswd,tdhrr,tdbsr)
    use params
    use partitioning
    use mpi
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer,intent(in) :: ichunk(3) ! portion of array to operate on 
    integer,intent(in) :: mu ! -3 <= mu <= 3
    ! Temp Border Partition Cube
    type(localpart),intent(in) :: tbpc(-1:1,-1:1,-1:1)
    ! Temp DSlash Split Work Done
    logical, intent(inout) :: tdsswd(-3:3,-1:1,-1:1,-1:1)
    ! Temp Dirac Halo Recv Requests
    integer,intent(inout) :: tdhrr(54)
    ! Temp Dirac Border Send Requests
    integer,intent(inout) :: tdbsr(54)


    integer :: chunk(2,3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: inn
    integer :: ierr


    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not.any(tdsswd(:,ichunk(1),ichunk(2),ichunk(3)))
    
    if(halo_to_wait_for.ne.0) then
      call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
    endif
    if(mu.eq.0)then
      call dslash_split_local(Phi,R,am,imass,chunk,init)
    else 
     if(mu.gt.0) then
        call dslash_split_nonlocal(Phi,R,u,chunk,mu,1,init)
      else if(mu.lt.0) then
        call dslash_split_nonlocal(Phi,R,u,chunk,-mu,-1,init)

      endif
    endif

    ! flagging work done
    tdsswd(mu,ichunk(1),ichunk(2),ichunk(3)) = .true.
    ! checking whether to send the partition already or not
    if(all(tdsswd(:,ichunk(1),ichunk(2),ichunk(3))))then
      tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
      do inn=1,tpart%nn
        ! clearing send requests
        call MPI_Wait(tdbsr(tpart%ahpss(inn)),MPI_STATUS_IGNORE,ierr)
        ! restarting send request
        call MPI_Start(tdbsr(tpart%ahpss(inn)),ierr)
      enddo
    endif

  end subroutine

  !pure subroutine dslash_split_nonlocal(Phi,R,u,chunk,mu,v,init)
  subroutine dslash_split_nonlocal(Phi,R,u,chunk,mu,v,init)
    use params
    use dirac
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on
    integer, intent(in) :: mu ! 1 <= mu <= 3
    integer, intent(in) :: v
    logical, intent(in) :: init
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

    if(init)then
      if(v.eq.1) then
        do idirac=1,4
          igork=gamin(mu,idirac)
          do it = td,tu
            do iy = yd,yu
              do ix = xd,xu
                Phi(:,ix,iy,it,idirac)= &
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
                Phi(:,ix,iy,it,idirac)= &
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
    else ! not init 
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
    endif

  end subroutine dslash_split_nonlocal 

  pure subroutine dslash_split_local(Phi,R,am,imass,chunk,init)
    use params
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on
    logical, intent(in) :: init
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
    if(init)then
      Phi(:,xd:xu,yd:yu,td:tu,:) = diag * R(:,xd:xu,yd:yu,td:tu,:)
    else! not init
      Phi(:,xd:xu,yd:yu,td:tu,:) = Phi(:,xd:xu,yd:yu,td:tu,:) + &
       & diag * R(:,xd:xu,yd:yu,td:tu,:)
    endif
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

! DSLASHD
  subroutine dslashd_split(Phi,R,u,am,imass,ichunk,mu,tbpc,tdsswd,tdhrr,tdbsr)
    use params
    use partitioning
    use mpi
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer,intent(in) :: ichunk(3) ! portion of array to operate on 
    integer,intent(in) :: mu ! -3 <= mu <= 3
    ! Temp Border Partition Cube
    type(localpart),intent(in) :: tbpc(-1:1,-1:1,-1:1)
    ! Temp DSlash Split Work Done
    logical, intent(inout) :: tdsswd(-3:3,-1:1,-1:1,-1:1)
    ! Temp Dirac Halo Recv Requests
    integer,intent(inout) :: tdhrr(54)
    ! Temp Dirac Border Send Requests
    integer,intent(inout) :: tdbsr(54)



    integer :: chunk(2,3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: inn
    integer :: ierr


    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not.any(tdsswd(:,ichunk(1),ichunk(2),ichunk(3)))
    
    if(halo_to_wait_for.ne.0) then
      call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
    endif
    if(mu.eq.0)then
      call dslashd_split_local(Phi,R,am,imass,chunk,init)
    else 
     if(mu.gt.0) then
        call dslashd_split_nonlocal(Phi,R,u,chunk,mu,1,init)
      else if(mu.lt.0) then
        call dslashd_split_nonlocal(Phi,R,u,chunk,-mu,-1,init)

      endif
    endif

    ! flagging work done
    tdsswd(mu,ichunk(1),ichunk(2),ichunk(3)) = .true. 
    ! checking whether to send the partition already or not
    if(all(tdsswd(:,ichunk(1),ichunk(2),ichunk(3))))then
      tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
      do inn=1,tpart%nn
        ! clearing send requests
        call MPI_Wait(tdbsr(tpart%ahpss(inn)),MPI_STATUS_IGNORE,ierr)
        ! restarting send request
        call MPI_Start(tdbsr(tpart%ahpss(inn)),ierr)
      enddo
    endif


  end subroutine

  pure subroutine dslashd_split_nonlocal(Phi,R,u,chunk,mu,v,init)
    use params
    use dirac
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on
    integer, intent(in) :: mu ! 1 <= mu <= 3
    integer, intent(in) :: v
    logical, intent(in) :: init
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

    if(init)then
      if(v.eq.1) then
        do idirac=1,4
          igork=gamin(mu,idirac)
          do it = td,tu
            do iy = yd,yu
              do ix = xd,xu
                Phi(:,ix,iy,it,idirac)= &
                  ! Wilson term (hermitian)
                &    -akappa*(u(ix,iy,it,mu) &
                  &              * R(:, ix+ixup, iy+iyup, it+itup, idirac)) &
                  ! Dirac term (antihermitian)
                &     - gamval(mu,idirac) * &
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
                Phi(:,ix,iy,it,idirac)= &
                  ! Wilson term (hermitian)
                &    -akappa*( conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &              * R(:, ix-ixup, iy-iyup, it-itup, idirac)) &
                  ! Dirac term (antihermitian)
                &     - gamval(mu,idirac) * &
                  &       (- conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
              enddo
            enddo
          enddo
        enddo
      endif
    else ! not init 
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
                &     - gamval(mu,idirac) * &
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
                &     - gamval(mu,idirac) * &
                  &       (- conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                  &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  end subroutine dslashd_split_nonlocal 

  pure subroutine dslashd_split_local(Phi,R,am,imass,chunk,init)
    use params
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on
    logical, intent(in) :: init
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
    if(init)then
      Phi(:,xd:xu,yd:yu,td:tu,:) = diag * R(:,xd:xu,yd:yu,td:tu,:)
    else ! not init
      Phi(:,xd:xu,yd:yu,td:tu,:) = Phi(:,xd:xu,yd:yu,td:tu,:) +&
        & diag * R(:,xd:xu,yd:yu,td:tu,:)
    endif

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


  ! A guess at the best ordering for computing pieces for work to do.
  subroutine get_dslash_work_ordering(tdswo,bbf)
    use params
    use comms
    implicit none
    ! Temp DSlash Work Ordering
    integer, intent(out) :: tdswo(4,27*7)
    ! Bulk Before Flag
    logical, intent(in) :: bbf
    ! Work Partition Count 
    integer :: wpc 
    integer :: ipx,ipy,ipt
    integer :: ips(3)
    integer :: ip2sum
    integer :: mu,musign,muabs
    logical :: needs_comms_before(-3:3,-1:1,-1:1,-1:1)
    ! Process Parity
    ! It is actually not defined for odd-sized grids, but there's nothing 
    ! better.
    integer :: pp

    ! ARBITRARY, the order could be done differently
    ! e.g., all chunks that do not need communication first and then 
    ! the ones that do.
    do ipt=-1,1
      do ipy=-1,1
        do ipx=-1,1
          ips = (/ipx,ipy,ipt/)
          ! determining first wheter or not a workload needs communications
          do mu=-3,3
            if(mu.eq.0) then
              needs_comms_before(mu,ipx,ipy,ipt) = .false.
            else
              musign = sign(1,mu)
              muabs = abs(mu)
              needs_comms_before(mu,ipx,ipy,ipt) = musign.eq.ips(muabs)
            endif
          enddo
        enddo
      enddo
    enddo


    wpc = 0

    if(bbf)then ! Taking care of the bulk first
      do mu=-3,3
        wpc = wpc + 1
        tdswo(:,wpc) = (/0,0,0,mu/)
      enddo
    endif

    ! selecting first the directions that don't need communications before
    ! vertices first, then edges, then faces (bulk treated separately)
    do ip2sum=3,1,-1
      do ipt=-1,1
        do ipy=-1,1
          do ipx=-1,1
            ! selecting vertices/edges/faces/bulk according to ip2sum
            if((ipx**2+ipy**2+ipt**2).eq.ip2sum)then
              do mu=-3,3
                ! selecting first the directions that don't need communications
                ! before
                if(.not.needs_comms_before(mu,ipx,ipy,ipt))then
                  wpc = wpc + 1
                  tdswo(:,wpc) = (/ipx,ipy,ipt,mu/)
                endif
              enddo
            endif
          enddo
        enddo
      enddo
    enddo

    ! and then the directions that do need communications before
    ! vertices first, then edges, then faces (bulk does not need communications
    ! before)
    ! trying to switch order in a red-black-grid-like fashion
    pp = 2*mod(ip_x+ip_y+ip_t,2)-1
    do ip2sum=3,1,-1
      do ipt=-pp,pp,pp
        do ipy=-pp,pp,pp
          do ipx=-pp,pp,pp
            ! selecting vertices/edges/faces/bulk according to ip2sum
            if((ipx**2+ipy**2+ipt**2).eq.ip2sum)then
              do mu=-3,3
                ! and then the directions that do need communications before
                if(needs_comms_before(mu,ipx,ipy,ipt))then
                  wpc = wpc + 1
                  tdswo(:,wpc) = (/ipx,ipy,ipt,mu/)
                endif
              enddo
            endif
          enddo
        enddo
      enddo
    enddo

  if(.not.bbf)then ! Taking care of the bulk last
    do mu=-3,3
      wpc = wpc + 1
      tdswo(:,wpc) = (/0,0,0,mu/)
    enddo
  endif



  end subroutine


end module dirac_split




 
