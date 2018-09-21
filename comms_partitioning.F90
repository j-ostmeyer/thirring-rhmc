module comms_partitioning
  use partitioning
  use mpi_f08

  type(MPI_Datatype) :: dirac_halo_dts(54)   ! Data TypeS
  type(MPI_Datatype) :: dirac_border_dts(26) ! Data TypeS

contains

  subroutine init_dirac_halo_types(tdhdts,thpl)
    use partitioning
    use mpi_f08
    type(MPI_Datatype),intent(out) :: tdhdts(54) ! Temp Dirac Halo Data TypeS
    type(localpart),intent(in)     :: thpl(54)   ! Temp Halo Partition List
    integer :: ih
    integer :: chunk(2,3)

    integer, dimension(5) :: sizes, subsizes, starts

    integer :: ierr

    sizes=(/kthird,ksizex_l+2,ksizey_l+2,ksizet_l+2,4/)
    subsizes=(/kthird,-1,-1,-1,4/)
    starts=(/0,-1,-1,-1,0/) !
    do ih=1,54
      chunk = thpl(ih)%chunk
      subsizes(2:4) = chunk(2,:) - chunk(1,:)
      ! For MPI_Type_create_subarray, indices start at 0 (because it's a C function)
      starts(2:4) = chunk(1,:)-1 
      call MPI_Type_create_subarray(5,sizes,subsizes,start,MPI_Order_Fortran,&
      &  MPI_Double_Complex, tdhdts(ih),ierr)
      call MPI_Type_Commit(tdhdts(ih),ierr)
    enddo
  end subroutine

  subroutine init_dirac_border_types(tdbdts,tbpl)
    use partitioning
    use mpi_f08
    type(MPI_Datatype),intent(out) :: tdbdts(26) ! Temp Dirac Border Data TypeS
    type(localpart),intent(in)     :: tbpl(26)   ! Temp Border Partition List
    integer :: ib
    integer :: chunk(2,3)

    integer, dimension(5) :: sizes, subsizes, starts

    integer :: ierr

    sizes=(/kthird,ksizex_l+2,ksizey_l+2,ksizet_l+2,4/)
    subsizes=(/kthird,-1,-1,-1,4/)
    starts=(/0,-1,-1,-1,0/) !
    do ib=1,54
      chunk = tbpl(ih)%chunk
      subsizes(2:4) = chunk(2,:) - chunk(1,:)
      ! For MPI_Type_create_subarray, indices start at 0 (because it's a C function)
      starts(2:4) = chunk(1,:)-1 
      call MPI_Type_create_subarray(5,sizes,subsizes,start,MPI_Order_Fortran,&
      &  MPI_Double_Complex, tdhdts(ih))
      call MPI_Type_Commit(tdbdts(ih),ierr)
    enddo
  end subroutine

! TODO: create persistent comm requests
! TODO: test!

end module comms_partitioning
