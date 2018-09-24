module dirac_split
! Versions of the dirac operator subroutines which operate only on a subset 
! of the lattice. 
  use params
  use comms
  implicit none
  save

  logical :: dslash_swd(-3:3,-1:1,-1:1,-1:1)  ! DSLASH Split Work Done
  logical :: dslashd_swd(-3:3,-1:1,-1:1,-1:1) ! DSLASHD Split Work Done

contains 

  subroutine check_work_and_start_transfer(tdswd,tbpc,ichunk,sreqs)
    use mpi_f08
    use partitioning
    implicit none
    !Temp DSLASH(D) Split Work Done
    logical, intent(in) :: tdswd(-3:3,-1:1,-1:1,-1:1) 
    ! Temp Border Partition Cube
    type(localpart), intent(in) :: tbpc(-1:1,-1:1,-1:1)
    integer,intent(in) :: ichunk(3) ! portion of array to operate on 
    type(MPI_Request) :: sreqs(54) ! persistent Send REQuestS

    integer :: nn,inn
    type(localpart) :: tpart
    integer :: ierr

    if(all(tdswd(:,ichunk(1),ichunk(2),ichunk(3))))then
      tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
      nn = tpart%nn
      do inn=1,nn
        call MPI_Start(sreqs(tpart%ahpss(inn)),ierr)
      enddo
    endif
  end subroutine

  subroutine dslash_split(Phi,R,u,am,imass,ichunk,mu,tbpc,tdsswd,tdhrr)
    use partitioning
    use mpi_f08
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
    type(MPI_Request) :: tdhrr(54)

    integer :: chunk(2,3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: ierr


    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not.any(tdsswd(:,ichunk(1),ichunk(2),ichunk(3)))

    if(mu.eq.0)then
      call dslash_split_local(Phi,R,am,imass,chunk,init)
    else 
      if(halo_to_wait_for.ne.0) then
        call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
      endif
      if(mu.gt.0) then
        call dslash_split_nonlocal(Phi,R,u,chunk,mu,2,init)
      else if(mu.lt.0) then
        call dslash_split_nonlocal(Phi,R,u,chunk,-mu,1,init)

      endif
    endif

    ! flagging work done
    tdsswd(mu,ichunk(1),ichunk(2),ichunk(3)) = .true.

  end subroutine

  pure subroutine dslash_split_nonlocal(Phi,R,u,chunk,mu,v,init)
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

  subroutine dslashd_split(Phi,R,u,am,imass,ichunk,mu,tbpc,tdsswd,tdhrr)
    use partitioning
    use mpi_f08
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
    type(MPI_Request) :: tdhrr(54)

    integer :: chunk(2,3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: ierr


    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not.any(tdsswd(:,ichunk(1),ichunk(2),ichunk(3)))

    if(mu.eq.0)then
      call dslashd_split_local(Phi,R,am,imass,chunk,init)
    else 
      if(halo_to_wait_for.ne.0) then
        call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
      endif
      if(mu.gt.0) then
        call dslashd_split_nonlocal(Phi,R,u,chunk,mu,2,init)
      else if(mu.lt.0) then
        call dslashd_split_nonlocal(Phi,R,u,chunk,-mu,1,init)

      endif
    endif

    ! flagging work done
    tdsswd(mu,ichunk(1),ichunk(2),ichunk(3)) = .true. 

  end subroutine


  pure subroutine dslashd_split_nonlocal(Phi,R,u,chunk,mu,v,init)
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

  end subroutine dslashd_split_nonlocal 

  pure subroutine dslash_split_local(Phi,R,am,imass,chunk,init)
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

  subroutine dslashd_split(Phi,R,u,am,imass,ichunk,mu,tbpc,tdsswd,tdhrr)
    use partitioning
    use mpi_f08
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
    type(MPI_Request) :: tdhrr(54)

    integer :: chunk(2,3)
    logical :: init
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: ierr


    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk
    halo_to_wait_for = tpart%ahpsr(mu)
    ! checking if some work on the partition has already been done
    init = .not.any(tdsswd(:,ichunk(1),ichunk(2),ichunk(3)))

    if(mu.eq.0)then
      call dslashd_split_local(Phi,R,am,imass,chunk,init)
    else 
      if(halo_to_wait_for.ne.0) then
        call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
      endif
      if(mu.gt.0) then
        call dslashd_split_nonlocal(Phi,R,u,chunk,mu,2,init)
      else if(mu.lt.0) then
        call dslashd_split_nonlocal(Phi,R,u,chunk,-mu,1,init)

      endif
    endif

    ! flagging work done
    tdsswd(mu,ichunk(1),ichunk(2),ichunk(3)) = .true. 

  end subroutine


  pure subroutine dslashd_split_nonlocal(Phi,R,u,chunk,mu,v,init)
    use dirac
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on
    integer, intent(in) :: mu
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
                phi(:,ix,iy,it,idirac)= &
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
                phi(:,ix,iy,it,idirac)= &
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
    else ! not init
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
    endif


  end subroutine dslashd_split_nonlocal

  pure subroutine dslashd_split_local(Phi,R,am,imass,chunk,init)
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

end module dirac_split





