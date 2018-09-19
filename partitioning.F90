module partitioning
  use params
  implicit none
  integer :: all_partitions(2,3,-2:2,-2:2,-2:2) ! with halo 
  integer :: all_local_partitions_neighbors(0:3,2,-1:1,-1:1,-1:1)
  ! [ no_neighbors, neigh1 , neigh2 , neigh3  ]
  ! [ ////////////, mpitag1, mpitag2, mpitag3 ] ! tags go from 1 to 9

  integer :: all_halo_partitions_neighbors(2,-2:2,-2:2,-2:2)
  ! [ neigh, tag ]

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

  subroutine get_all_local_partitions_neighbors(temp_all_local_partitions_neighbors)
    use comms
    use mpi_f08
    integer,intent(out) :: temp_all_local_partitions_neighbors(0:3,2,-1:1,-1:1,-1:1)
    integer :: ipx,ipy,ipt
    integer :: ipn,ipw ! ipw not actually used
    integer :: neighbour_count
    integer :: ipxs(3)
    integer :: neighidx, idir
    integer :: facecoord(2),tag
    integer :: ierr

    do ipx=-1,1
      do ipy=-1,1
        do ipt=-1,1
          ipxs = (/ ipx, ipy, ipt /)
          neighbour_count = sum(ipxs**2)

          temp_all_local_partitions_neighbors(0,1,ipx,ipy,ipt) = neighbour_count
          temp_all_local_partitions_neighbors(0,2,ipx,ipy,ipt) = -1 ! unused
          neighidx = 1
          do idir=1,3
            if(ipxs(idir).ne.0) then
              ! getting neighbour
              call MPI_Cart_Shift(comm, 0, ipxs(idir), ipw, ipn,ierr)
              temp_all_local_partitions_neighbors(neighidx,1,ipx,ipy,ipt) = ipn

              ! computing tag
              facecoord(1:idir-1)=ipxs(1:idir-1)
              facecoord(idir:2)=ipxs(idir+1:3)
              tag = 1+(facecoord(1)+1)+(facecoord(2)+1)*3
              
              temp_all_local_partitions_neighbors(neighidx,2,ipx,ipy,ipt) = tag

              neighidx = neighidx + 1 
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine

  subroutine get_all_halo_partitions_neighbors(temp_all_halo_partitions_neighbors)
    use comms
    use mpi_f08
    integer, intent(out) :: temp_all_halo_partitions_neighbors(2,-2:2,-2:2,-2:2)
    integer :: ipx,ipy,ipt
    integer :: ipn,ipw ! ipw not actually used
    integer :: ipxs(3)
    integer :: idir
    integer :: facecoord(2),tag
    integer :: ierr
    integer :: discr


    temp_all_halo_partitions_neighbors = -1

    do ipx=-2,2
      do ipy=-2,2
        do ipt=-2,2
          ipxs = (/ ipx, ipy, ipt /)
          discr = sum(ipxs**2)
          if((discr.ge.4).and.(discr.le.6)) then 
            ! then we are in the halo. but not on the edges or on the vertices
            do idir=1,3
              if(ipxs(idir)**2.eq.4) then
                ! getting neighbour
                call MPI_Cart_Shift(comm, 0, ipxs(idir)/2, ipw, ipn,ierr)
                temp_all_halo_partitions_neighbors(1,ipx,ipy,ipt) = ipn

                ! computing tag
                facecoord(1:idir-1)=ipxs(1:idir-1)
                facecoord(idir:2)=ipxs(idir+1:3)
                tag = 1+(facecoord(1)+1)+(facecoord(2)+1)*3

                temp_all_halo_partitions_neighbors(2,ipx,ipy,ipt) = tag

              endif
            enddo
          endif
        enddo
      enddo
    enddo

  end subroutine


  subroutine init_partitions_and_neighs()
    call get_all_partitions(all_partitions)
    call get_all_local_partitions_neighbors(all_local_partitions_neighbors)
    call get_all_halo_partitions_neighbors(all_halo_partitions_neighbors)
  end subroutine


end module partitioning
