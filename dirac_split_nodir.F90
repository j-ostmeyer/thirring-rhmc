module dirac_split_nodir
! Versions of the dirac operator subroutines which operate only on a subset 
! of the lattice. 
  use params
  use comms
  implicit none
  save

  !logical :: dslash_swd(-1:1,-1:1,-1:1)  ! DSLASH Split Work Done
  !logical :: dslashd_swd(-1:1,-1:1,-1:1) ! DSLASHD Split Work Done

  ! list of partitions, indexed by ipx,ipy and ipt partition coordinates, 
  ! and mu.
  ! with  -1 <= ip[xyt] <= 1, -3 <= mu <= 3
  ! in principle it can also depend on the subroutine they are used in 
  ! (qmrherm or congrad)
  integer :: dslash_work_ordering(3,27)
  integer :: dslashd_work_ordering(3,27)

contains 

! DSLASH 

  subroutine dslash_split(Phi,R,u,am,imass,ichunk,tbpc,tdhrr,tdbsr)
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
    ! Temp Border Partition Cube
    type(localpart),intent(in) :: tbpc(-1:1,-1:1,-1:1)
    ! Temp Dirac Halo Recv Requests
    integer,intent(inout) :: tdhrr(54)
    ! Temp Dirac Border Send Requests
    integer,intent(inout) :: tdbsr(54)


    integer :: chunk(2,3)
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: inn
    integer :: ierr

    integer :: mu


    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk

    do mu=-3,3
      halo_to_wait_for = tpart%ahpsr(mu)
      if(halo_to_wait_for.ne.0) then
        call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
      endif
    enddo

    call dslash_split_work(Phi,R,am,imass,u,chunk)

    ! send the partition
    do inn=1,tpart%nn
      ! clearing send requests
      call MPI_Wait(tdbsr(tpart%ahpss(inn)),MPI_STATUS_IGNORE,ierr)
      ! restarting send request
      call MPI_Start(tdbsr(tpart%ahpss(inn)),ierr)
    enddo

  end subroutine

  subroutine dslash_split_work(Phi,R,am,imass,u,chunk)
    use params
    use dirac
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: am
    integer, intent(in) :: imass
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on

    integer :: mu
    integer :: xd,xu,yd,yu,td,tu ! portion of array to operate on
    integer :: ixup, iyup, itup, ix, iy, it, idirac, igork
    real(dp) :: diag
    complex(dp) :: zkappa

    xd=chunk(1,1)
    xu=chunk(2,1)
    yd=chunk(1,2)
    yu=chunk(2,2)
    td=chunk(1,3)
    tu=chunk(2,3)

    diag=(3.0-am3)+1.0
    do it = td,tu
      do iy = yd,yu
        do ix = xd,xu
          ! local term (initialization)
          do idirac=1,4
            Phi(:,ix,iy,it,idirac) = diag * R(:,ix,iy,it,idirac)

            ! nonlocal term
            do mu=1,3
              ixup = kdelta(1, mu)
              iyup = kdelta(2, mu)
              itup = kdelta(3, mu)
              igork=gamin(mu,idirac)

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
           
            enddo !do mu=1,3
          enddo !do idirac=1,4

          ! Local term
          !      
          !  s-like term exploiting projection
          Phi(1:kthird-1, ix, iy, it, 3:4) &
            & = Phi(1:kthird-1, ix, iy, it, 3:4) &
            & - R(2:kthird, ix, iy, it, 3:4)
          Phi(2:kthird, ix, iy, it, 1:2) &
            & = Phi(2:kthird, ix, iy, it, 1:2) &
            & - R(1:kthird-1, ix, iy, it, 1:2)
          !
          !  Mass term (couples the two walls unless imass=5)
          if (imass.eq.1) then
            zkappa=cmplx(am,0.0)
            Phi(kthird, ix, iy, it, 3:4) = &
              & Phi(kthird, ix, iy, it, 3:4) &
              & + zkappa * R(1, ix, iy, it, 3:4)
            Phi(1, ix, iy, it, 1:2) = &
              & Phi(1, ix, iy, it, 1:2) + &
              & zkappa * R(kthird, ix, iy, it, 1:2)
          elseif (imass.eq.3) then
            zkappa=cmplx(0.0,-am)
            Phi(kthird, ix, iy, it, 3:4) = &
              & Phi(kthird, ix, iy, it, 3:4) &
              & - zkappa * R(1, ix, iy, it, 3:4)
            Phi(1, ix, iy, it, 1:2) = &
              & Phi(1, ix, iy, it, 1:2) &
              & + zkappa * R(kthird, ix, iy, it, 1:2)
          elseif (imass.eq.5) then
            zkappa=cmplx(0.0,-am)
            Phi(kthird, ix, iy, it, 3:4) = &
              & Phi(kthird, ix, iy, it, 3:4) &
              & - zkappa * R(kthird, ix, iy, it, 1:2)
            Phi(1, ix, iy, it, 1:2) = &
              & Phi(1, ix, iy, it, 1:2) &
              & - zkappa * R(1, ix, iy, it, 3:4)
          endif
 
        enddo! do ix = xd,xu
      enddo! do iy = yd,yu
    enddo! do it = td,tu
 
  end subroutine


  subroutine dslashd_split(Phi,R,u,am,imass,ichunk,tbpc,tdhrr,tdbsr)
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
    ! Temp Border Partition Cube
    type(localpart),intent(in) :: tbpc(-1:1,-1:1,-1:1)
    ! Temp Dirac Halo Recv Requests
    integer,intent(inout) :: tdhrr(54)
    ! Temp Dirac Border Send Requests
    integer,intent(inout) :: tdbsr(54)


    integer :: chunk(2,3)
    integer :: halo_to_wait_for
    type(localpart) :: tpart
    integer :: inn
    integer :: ierr

    integer :: mu

    tpart = tbpc(ichunk(1),ichunk(2),ichunk(3))
    chunk = tpart%chunk
    do mu=-3,3
      halo_to_wait_for = tpart%ahpsr(mu)
      if(halo_to_wait_for.ne.0) then
        call MPI_Wait(tdhrr(halo_to_wait_for),MPI_STATUS_IGNORE,ierr)
      endif
    enddo

    call dslashd_split_work(Phi,R,am,imass,u,chunk)
    ! send the partition
    do inn=1,tpart%nn
      ! clearing send requests
      call MPI_Wait(tdbsr(tpart%ahpss(inn)),MPI_STATUS_IGNORE,ierr)
      ! restarting send request
      call MPI_Start(tdbsr(tpart%ahpss(inn)),ierr)
    enddo

  end subroutine

  subroutine dslashd_split_work(Phi,R,am,imass,u,chunk)
    use params
    use dirac
    implicit none
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: am
    integer, intent(in) :: imass
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    integer, intent(in) :: chunk(2,3) ! portion of array to operate on

    integer :: mu
    integer :: xd,xu,yd,yu,td,tu ! portion of array to operate on
    integer :: ixup, iyup, itup, ix, iy, it, idirac, igork
    real(dp) :: diag
    complex(dp) :: zkappa

    xd=chunk(1,1)
    xu=chunk(2,1)
    yd=chunk(1,2)
    yu=chunk(2,2)
    td=chunk(1,3)
    tu=chunk(2,3)


    diag=(3.0-am3)+1.0
    do it = td,tu
      do iy = yd,yu
        do ix = xd,xu
          do idirac=1,4
            ! local term (initialisation)
            Phi(:,ix,iy,it,idirac) = diag * R(:,ix,iy,it,idirac)
            ! nonlocal term
            do mu=1,3
              igork=gamin(mu,idirac)
              ixup = kdelta(1, mu)
              iyup = kdelta(2, mu)
              itup = kdelta(3, mu)

              Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac) &
                !   wilson term (hermitian)
                &    - akappa * (u(ix,iy,it,mu) &
                &              * R(:, ix+ixup, iy+iyup, it+itup, idirac) &
                &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &              * R(:, ix-ixup, iy-iyup, it-itup, idirac)) &
                !   dirac term (antihermitian)
                &    - gamval(mu,idirac) * &
                &       (u(ix,iy,it,mu) &
                &         * R(:, ix+ixup, iy+iyup, it+itup, igork) &
                &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
            enddo! do mu=1,3
          enddo! do idirac=1,4
 
          ! local term
          !   s-like term exploiting projection
          Phi(1:kthird-1, ix, iy, it, 1:2) &
            & = Phi(1:kthird-1, ix, iy, it, 1:2) &
            & - R(2:kthird, ix, iy, it, 1:2)
          Phi(2:kthird, ix, iy, it, 3:4) &
            & = Phi(2:kthird, ix, iy, it, 3:4) &
            & - R(1:kthird-1, ix, iy, it, 3:4)
          !
          !   Mass term (couples the two walls unless imass=5) 
          if(imass.eq.1)then
            zkappa=cmplx(am,0.0)
            Phi(kthird, ix, iy, it, 1:2) = &
              & Phi(kthird, ix, iy, it, 1:2) &
              & + zkappa * R(1, ix, iy, it, 1:2)
            Phi(1, ix, iy, it, 3:4) = &
              & Phi(1, ix, iy, it, 3:4) &
              & + zkappa * R(kthird, ix, iy, it, 3:4)
          elseif(imass.eq.3)then
            zkappa = cmplx(0.0,am)
            Phi(kthird, ix, iy, it, 1:2) = &
              & Phi(kthird, ix, iy, it, 1:2) &
              & + zkappa * R(1, ix, iy, it, 1:2)
            Phi(1, ix, iy, it, 3:4) = &
              & Phi(1, ix, iy, it, 3:4) &
              & - zkappa * R(kthird, ix, iy, it, 3:4)
          elseif(imass.eq.5)then
            zkappa = cmplx(0.0,am)
            Phi(kthird, ix, iy, it, 1:2) = &
              & Phi(kthird, ix, iy, it, 1:2) &
              & - zkappa * R(kthird, ix, iy, it, 3:4)
            Phi(1, ix, iy, it, 3:4) = &
              & Phi(1,ix, iy, it, 3:4) &
              & - zkappa * R(1, ix, iy, it, 1:2)
          endif
        enddo
      enddo
    enddo
  end subroutine

  ! A guess at the best ordering for computing pieces for work to do.
  subroutine get_dslash_work_ordering(tdswo,bbf)
    use params
    use comms
    implicit none
    ! Temp DSlash Work Ordering
    integer, intent(out) :: tdswo(3,27)
    ! Bulk Before Flag
    logical, intent(in) :: bbf
    ! Work Partition Count 
    integer :: wpc 
    integer :: ipx,ipy,ipt
    integer :: ip2sum
    logical :: needs_comms_before(-1:1,-1:1,-1:1)
    ! Process Parity
    ! It is actually not defined for odd-sized grids, but there's nothing 
    ! better.
    integer :: pp

    ! ARBITRARY, the order could be done differently
    ! e.g., all chunks that do not need communication first and then 
    ! the ones that do.
    needs_comms_before = .true.
    needs_comms_before(0,0,0) = .false.

    wpc = 0

    if(bbf)then ! Taking care of the bulk first
      wpc = wpc + 1
      tdswo(:,wpc) = (/0,0,0/)
    endif

    ! and then the surfaces, that do need communications before
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
              wpc = wpc + 1
              tdswo(:,wpc) = (/ipx,ipy,ipt/)
            endif
          enddo
        enddo
      enddo
    enddo

    if(.not.bbf)then ! Taking care of the bulk last
      wpc = wpc + 1
      tdswo(:,wpc) = (/0,0,0/)
    endif



  end subroutine


end module dirac_split_nodir




 
