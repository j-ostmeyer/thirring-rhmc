module comms_partitioning
  use partitioning
  use mpi
  implicit none

  integer :: dirac_halo_dts(54)   ! Data TypeS
  integer :: dirac_border_dts(26) ! Data TypeS

  !integer :: dirac_border_send_reqs(54)
  !integer :: dirac_halo_recv_reqs(54)

contains

  subroutine init_dirac_halo_types(tdhdts, thpl)
    use params
    use partitioning
    use mpi
    implicit none
    integer, intent(out) :: tdhdts(54) ! Temp Dirac Halo Data TypeS
    type(halopart), intent(in)     :: thpl(54)   ! Temp Halo Partition List
    integer :: ih
    integer :: chunk(2, 3)

    integer, dimension(5) :: sizes, subsizes, starts

    integer :: ierr

    sizes = (/kthird, ksizex_l + 2, ksizey_l + 2, ksizet_l + 2, 4/)
    subsizes = (/kthird, -1, -1, -1, 4/)
    starts = (/0, -1, -1, -1, 0/) !
    do ih = 1, 54
      chunk = thpl(ih)%chunk
      subsizes(2:4) = chunk(2, :) - chunk(1, :) + 1
      ! For MPI_Type_create_subarray, indices start at 0 (because it's a C function)
      starts(2:4) = chunk(1, :) ! arrays are indexed from 0 as chunks
      call MPI_Type_create_subarray(5, sizes, subsizes, starts, MPI_Order_Fortran,&
          &  MPI_Double_Complex, tdhdts(ih), ierr)
      call MPI_Type_Commit(tdhdts(ih), ierr)
    enddo
  end subroutine

  subroutine init_dirac_border_types(tdbdts, tbpl)
    use params
    use partitioning
    use mpi
    implicit none
    integer, intent(out) :: tdbdts(26) ! Temp Dirac Border Data TypeS
    type(localpart), intent(in)     :: tbpl(26)   ! Temp Border Partition List
    integer :: ib
    integer :: chunk(2, 3)

    integer, dimension(5) :: sizes, subsizes, starts

    integer :: ierr

    sizes = (/kthird, ksizex_l + 2, ksizey_l + 2, ksizet_l + 2, 4/)
    subsizes = (/kthird, -1, -1, -1, 4/)
    starts = (/0, -1, -1, -1, 0/) !
    do ib = 1, 26
      chunk = tbpl(ib)%chunk
      subsizes(2:4) = chunk(2, :) - chunk(1, :) + 1
      ! For MPI_Type_create_subarray, indices start at 0 (because it's a C function)
      starts(2:4) = chunk(1, :) ! arrays are indexed from 0 as chunks
      call MPI_Type_create_subarray(5, sizes, subsizes, starts, MPI_Order_Fortran,&
          &  MPI_Double_Complex, tdbdts(ib), ierr)
      call MPI_Type_Commit(tdbdts(ib), ierr)
    enddo
  end subroutine

  ! CREATE Dirac Persistent Send REQuestS
  subroutine create_dpsreqs(sreqs, bufts, tsbdts, tbpl)
    use mpi
    use comms
    use params
    use partitioning
    implicit none
    integer, intent(out) :: sreqs(54)! Send REQuestS
    !BUFfer To Send
    complex(dp), intent(in) :: bufts(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: tsbdts(26) ! Temp Send Border Data TypeS
    type(localpart), intent(in) :: tbpl(26) ! Temp Border Partition List
    integer :: ibp ! Index Border Partition
    integer :: ibhas ! Index Border Halo ASsociation

    integer :: ierr

    do ibp = 1, 26
      do ibhas = 1, tbpl(ibp)%nn
        call MPI_Send_init(bufts, 1, tsbdts(ibp),&
          & tbpl(ibp)%nns(ibhas), tbpl(ibp)%tags(ibhas),&
          & comm, sreqs(tbpl(ibp)%ahpss(ibhas)), ierr)
      enddo
    enddo
  end subroutine

  subroutine init_dirac_hb_types()
    use partitioning
    call init_dirac_halo_types(dirac_halo_dts, halo_partitions_list)
    call init_dirac_border_types(dirac_border_dts, border_partitions_list)
  end subroutine

  ! for convenience. GLOBAL data structures must be initialised first!
  subroutine get_dirac_sendreqs(sreqs, bufts)
    use params
    use partitioning
    implicit none
    integer, intent(out) :: sreqs(54)! Send REQuestS
    !BUFfer To Send
    complex(dp), intent(in) :: bufts(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call create_dpsreqs(sreqs, bufts, dirac_border_dts, border_partitions_list)
  end subroutine

  ! CREATE Dirac Persistent Rrecv REQuestS
  subroutine create_dprreqs(rreqs, buftr, trhdts, thpl)
    use mpi
    use comms
    use params
    use partitioning
    implicit none
    integer, intent(out) :: rreqs(54)! Recv REQuestS
    !BUFfer To Recv
    complex(dp), intent(in) :: buftr(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: trhdts(54) ! Temp Recv Halo Data TypeS
    type(halopart), intent(in) :: thpl(54) ! Temp Halo Partition List
    integer :: ihp ! Index Halo Partition

    integer :: ierr

    do ihp = 1, 54
      call MPI_Recv_init(buftr, 1, trhdts(ihp),&
        & thpl(ihp)%nn, thpl(ihp)%tag,&
        & comm, rreqs(ihp), ierr)
    enddo
  end subroutine

  ! for convenience. GLOBAL data structures must be initialised first!
  subroutine get_dirac_recvreqs(rreqs, buftr)
    use params
    use partitioning
    implicit none
    integer, intent(out) :: rreqs(54)! Recv REQuestS
    !BUFfer To Recv
    complex(dp), intent(in) :: buftr(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call create_dprreqs(rreqs, buftr, dirac_halo_dts, halo_partitions_list)
  end subroutine

end module comms_partitioning
