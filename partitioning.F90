module partitioning
  use params
  implicit none

  integer :: all_partitions(2,3,-2:2,-2:2,-2:2) ! with halo 

  type localpart
    integer :: chunk(2,3)
    integer :: nn ! number of neighbours the partition must be sent to (0..3)
    integer :: nns(3) ! ranks of the nearest neighbours
    integer :: tags(3) ! tags for the partition in all sends
  end type localpart

  type halopart
    integer :: chunk(2,3)
    integer :: nn ! rank of the neares neighbour
    integer :: tag ! tag to the partition
  end type halopart

  ! partitions as lists (contains only useful ones, but difficult to access)
  type(localpart) border_partitions_list(26)
  type(halopart) halo_partitions_list(54)

  ! partitions as a cube (contains also useless ones, but easier to access)
  type(localpart) border_partitions_cube(-1:1,-1:1,-1:1)
  type(halopart) halo_partitions_cube(-2:2,-2:2,-2:2)

  !conversions cube->list
  integer :: border_cl(-1:1,-1:1,-1:1)
  integer :: halo_cl(-2:2,-2:2,-2:2)
  !conversions list->cube
  integer :: border_lc(3,26)
  integer :: halo_lc(3,54)

  type bhas ! Border partition - Halo partition ASsociation (dictionary entry)
    ! Number of Halo Partition associated with the current border partition
    integer :: nhp 
    ! index (in halo_partitions_list) of the associated Halo Partitions
    integer :: hps(3) 
  end type bhas

  type(bhas) bhass(26) ! Border partition -> Halo partitions ASsociationS
 
contains
  ! Partition: a (2,3) array of integers with the contents
  ! ((xmin,xmax),(ymin,ymax),(tmin,tmax))
  ! expands first partition to include the second one. Just translates the 
  ! lower and upper limits.

  subroutine expand_partition(p1,p2)
    integer, intent(in) :: p2(2,3)
    integer, intent(inout) :: p1(2,3)
    integer :: idir

    do idir=1,3
      p1(1,idir) = min(p1(1,idir),p2(1,idir))
      p1(2,idir) = max(p1(2,idir),p2(2,idir))
    enddo

  end subroutine

  ! create an array with all the start-end indices for all partitions in
  ! all dimensions. 
  ! the partitions 
  subroutine get_border3geometry(border3geometry)
    integer,intent(out) :: border3geometry(2,-2:2,3)
    integer :: locdim(3), idimension
    integer :: portion_starts(-2:3),iportion

    locdim(1) = ksizex_l
    locdim(2) = ksizey_l
    locdim(3) = ksizet_l


    do idimension=1,3
      portion_starts(-2) = 0                   ! start halo
      portion_starts(-1) = 1                   ! start border
      portion_starts(0) = 2                    ! startbulk
      portion_starts(1) = locdim(idimension)   ! start border
      portion_starts(2) = locdim(idimension)+1 ! start halo
      portion_starts(3) = locdim(idimension)+2 ! beyond
      do iportion = -2,2
        border3geometry(1,iportion,idimension) = portion_starts(iportion)
        border3geometry(2,iportion,idimension) = portion_starts(iportion+1)-1
      enddo
    enddo

  end subroutine

  ! create a set of partitions
  subroutine get_all_partitions(temp_all_partitions)
    integer, intent(out) :: temp_all_partitions(2,3,-2:2,-2:2,-2:2)
    integer :: border3geometry(2,-2:2,3)
    integer :: ipx,ipy,ipt,iside
    integer :: ips(3)

    call get_border3geometry(border3geometry)
    do ipt=-2,2
      do ipy=-2,2
        do ipx=-2,2
          ips = (/ ipx,ipy,ipt /)
          do iside=1,3
            temp_all_partitions(:,iside,ipx,ipy,ipt) = & 
              & border3geometry(:,ips(iside),iside)
          enddo
        enddo
      enddo
    enddo

  end subroutine

  subroutine get_all_local_partitions_neighbors(tbp_cube,tbp_list,tb_cl,tb_lc)
    use comms
    use mpi_f08
    type(localpart),intent(out) :: tbp_cube(-1:1,-1:1,-1:1)!Temp Border Partitions
    type(localpart),intent(out) :: tbp_list(26)
    integer,intent(out) :: tb_cl(-1:1,-1:1,-1:1) ! Temp Border Cube->List
    integer,intent(out) :: tb_lc(3,26)           ! Temp Border List->Cube
    type(localpart) :: tpart
    integer :: tap(2,3,-2:2,-2:2,-2:2) ! temp_all_partitions
    integer :: ipx,ipy,ipt
    integer :: ipn,ipw ! ipw not actually used
    integer :: ips(3)
    integer :: neighidx, idir
    integer :: facecoord(2),tag
    integer :: ierr
    integer :: tpartcount

    tpartcount = 0
    tb_cl = 0
    tb_lc = -5555 ! an absurd value (everything shoulbe overwritten)
    call get_all_partitions(tap)
    do ipx=-1,1
      do ipy=-1,1
        do ipt=-1,1
          tpart%chunk = tap(:,:,ipx,ipy,ipt)
          ips = (/ ipx, ipy, ipt /)
 
          tpart%nn = sum(ips**2) ! neighbour count

          neighidx = 1
          do idir=1,3
            if(ips(idir).ne.0) then ! we have a neigbouring rank
              ! getting neighbour in direction idir
              call MPI_Cart_Shift(comm, idir-1, ips(idir), ipw, ipn,ierr)
              tpart%nns(neighidx) = ipn

              ! computing tag
              facecoord(1:idir-1)=ips(1:idir-1)
              facecoord(idir:2)=ips(idir+1:3)
              tag = 1+(facecoord(1)+1)+(facecoord(2)+1)*3
              
              tpart%tags(neighidx) = tag

              neighidx = neighidx + 1 
            endif
          enddo
          tbp_cube(ipx,ipy,ipt) = tpart
          if(tpart%nn.ne.0) then
            tpartcount = tpartcount+1
            tbp_list(tpartcount) = tpart
            tb_cl(ipx,ipy,ipt) = tpartcount 
            tb_lc(:,tpartcount) = (/ipx,ipy,ipt/)
          endif
        enddo
      enddo
    enddo
  end subroutine

  subroutine get_all_halo_partitions_neighbors(thp_cube,thp_list,th_cl,th_lc)
    use comms
    use mpi_f08
    type(halopart),intent(out) :: thp_cube(-2:2,-2:2,-2:2)!Temp Halo Partitions
    type(halopart),intent(out) :: thp_list(54)             
    integer, intent(out) :: th_cl(-2:2,-2:2,-2:2) ! Temp Halo Cube->List
    integer, intent(out) :: th_lc(3,54)           ! Temp Halo List->Cube
    type(halopart) :: tpart
    integer :: tap(2,3,-2:2,-2:2,-2:2) ! Temp All Partitions
    integer :: ipx,ipy,ipt
    integer :: ipn,ipw ! ipw not actually used
    integer :: ipxs(3)
    integer :: idir
    integer :: facecoord(2),tag
    integer :: ierr
    integer :: discr
    integer :: tpartcount

    tpartcount = 0
    th_cl = 0
    th_lc = -5555
    call get_all_partitions(tap)
    do ipx=-2,2
      do ipy=-2,2
        do ipt=-2,2
          ipxs = (/ ipx, ipy, ipt /)
          discr = sum(ipxs**2)
          if((discr.ge.4).and.(discr.le.6)) then 
            tpart%chunk = tap(:,:,ipx,ipy,ipt)
            ! then we are in the halo. but not on the edges or on the vertices
            do idir=1,3 ! scanning for ipdx(idir)=+-2
              if(ipxs(idir)**2.eq.4) then 
                ! getting neighbour
                call MPI_Cart_Shift(comm, idir-1, ipxs(idir)/2, ipw, ipn,ierr)
                tpart%nn = ipn

                ! computing tag
                facecoord(1:idir-1)=ipxs(1:idir-1)
                facecoord(idir:2)=ipxs(idir+1:3)
                tag = 1+(facecoord(1)+1)+(facecoord(2)+1)*3

                tpart%tag = tag
                write(7+ip_global,"(5I4)"),ipx,ipy,ipt,ipn,tag
              endif
            enddo
            thp_cube(ipx,ipy,ipt) = tpart
            tpartcount = tpartcount+1
            thp_list(tpartcount) = tpart

            th_cl(ipx,ipy,ipt) = tpartcount
            th_lc(:,tpartcount) = (/ipx,ipy,ipt/)
          endif
        enddo
      enddo
    enddo

  end subroutine

  subroutine get_border_halo_associations(tbhass,tb_lc,th_cl)
    type(bhas),intent(out) :: tbhass(26)        ! Temp bhass
    integer,intent(in) :: tb_lc(3,26)           ! Temp Border List->Cube
    integer,intent(in) :: th_cl(-2:2,-2:2,-2:2) ! Temp Halo Cube->List
    integer :: ips(3)
    integer :: hips(3) ! halo-ips
    integer :: idir
    integer :: ibp
    integer :: neighidx

    do ibp=1,26
      ips = tb_lc(:,ibp)
      tbhass(ibp)%nhp = sum(ips**2)
      neighidx = 0
      do idir=1,3
        if(ips(idir).ne.0) then ! we have a neigbouring rank
          hips = ips
          hips(idir) = -2*ips(idir) ! going on the halo on the opposite side
          neighidx = neighidx+1
          tbhass(ibp)%hps(neighidx) = th_cl(hips(1),hips(2),hips(3))
        endif
      enddo
    enddo
  end subroutine


  subroutine init_partitions_and_neighs()
    call get_all_partitions(all_partitions)
    call get_all_local_partitions_neighbors(border_partitions_cube,&
      &    border_partitions_list,border_cl,border_lc)
    call get_all_halo_partitions_neighbors(halo_partitions_cube,&
      &    halo_partitions_list,halo_cl,halo_lc)
    call get_border_halo_associations(bhass,border_lc,halo_cl)
  end subroutine


end module partitioning
