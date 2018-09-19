module locvol_partitioning
  use params
  implicit none
  integer :: united_border_partitions(2,3,2,3)
  integer :: bulk_partition(2,3)

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
    integer,intent(out) :: border3geometry(2,-1:1,3)
    integer :: locdim(3), idimension
    integer :: portion_starts(-1:2),iportion
    integer :: halo_size 

    locdim(1) = ksizex_l
    locdim(2) = ksizey_l
    locdim(3) = ksizet_l

    halo_size=1

    do idimension=1,3
      portion_starts(-1) = 1
      portion_starts(0) = halo_size+1
      portion_starts(1) = locdim(idimension)-halo_size+1
      portion_starts(2) = locdim(idimension)+1 ! already beyond
      do iportion = -1,1
        border3geometry(1,iportion,idimension) = portion_starts(iportion)
        border3geometry(2,iportion,idimension) = portion_starts(iportion+1)-1
      enddo
    enddo

  end subroutine

  ! create a set of partitions
  subroutine get_all_partitions(border_partitions)
    integer, intent(out) :: border_partitions(2,3,-1:1,-1:1,-1:1)
    integer :: border3geometry(2,-1:1,3)
    integer :: ipx,ipy,ipt,iside
    integer :: ips(3)

    call get_border3geometry(border3geometry)
    do ipt=-1,1
      do ipy=-1,1
        do ipx=-1,1
          ips = (/ ipx,ipy,ipt /)
          do iside=1,3
            border_partitions(:,iside,ipx,ipy,ipt) = & 
              & border3geometry(:,ips(iside),iside)
          enddo
        enddo
      enddo
    enddo

  end subroutine

  subroutine get_united_partitions(temp_united_border_partitions, &
    &               temp_bulk_partition)
    integer, intent(out) :: temp_united_border_partitions(2,3,2,3)
    integer, intent(out) :: temp_bulk_partition(2,3)
    integer ::  border_partitions(2,3,-1:1,-1:1,-1:1)
    integer :: ipx,ipy
    integer :: part(2,3)
    integer :: ud,udtgt

    call get_all_partitions(border_partitions)

    ! T border
    do ud=-1,1,2
      part = border_partitions(:,:,-1,-1,ud)
      do ipx=-1,1
        do ipy=-1,1
          call expand_partition(part,border_partitions(:,:,ipx,ipy,ud))
        enddo
      enddo
      udtgt = 1+(ud+1)/2
      temp_united_border_partitions(:,:,udtgt,3) = part
    enddo
    ! Y border
    do ud=-1,1,2
      part = border_partitions(:,:,-1,ud,0)
      do ipx=-1,1
        call expand_partition(part,border_partitions(:,:,ipx,ud,0))
      enddo
      udtgt = 1+(ud+1)/2
      temp_united_border_partitions(:,:,udtgt,2) = part
    enddo
    ! X xborder
    do ud=-1,1,2
      part = border_partitions(:,:,ud,0,0)
      udtgt = 1+(ud+1)/2
      temp_united_border_partitions(:,:,udtgt,1) = part
    enddo

    temp_bulk_partition = border_partitions(:,:,0,0,0)

  end subroutine

  subroutine init_partitions()
    call get_united_partitions(united_border_partitions,bulk_partition)
  end subroutine


end module locvol_partitioning
