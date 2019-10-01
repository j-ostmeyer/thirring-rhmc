module comms6
  use comms_common
  use params
  use mpi
  implicit none

  integer :: halo_6_xup_send(4:4, 1:25), halo_6_xdn_send(4:4, 1:25)
  integer :: halo_6_xup_recv(4:4, 1:25), halo_6_xdn_recv(4:4, 1:25)
  integer :: halo_6_yup_send(4:4, 1:25), halo_6_ydn_send(4:4, 1:25)
  integer :: halo_6_yup_recv(4:4, 1:25), halo_6_ydn_recv(4:4, 1:25)
  integer :: halo_6_tup_send(4:4, 1:25), halo_6_tdn_send(4:4, 1:25)
  integer :: halo_6_tup_recv(4:4, 1:25), halo_6_tdn_recv(4:4, 1:25)

contains
  !***********************************************************************
  !   Define MPI datatypes to allow halo transfer
  !***********************************************************************

  ! Initialise a single halo type.
  ! direction: direction to be communicated; x=0, y=1, t=2
  ! position: coordinate in direction of the region
  ! i5, i6: size in extra directions
  ! datatype: datatype to transfer
  ! typetarget: array element to put type into;
  !   should be halo_{4,4_real,5,6}_{x,y,t}{up,dn}_{send,recv}(size{4,5}[, size6])

  subroutine init_single_halo_type_6(direction, position, size5, size6, typetarget)
    integer, intent(in) :: direction, position, size5, size6
    integer, intent(out) :: typetarget
    integer, dimension(6) :: sizes, subsizes, starts
    integer :: ierr

    sizes = (/kthird, ksizex_l + 2, ksizey_l + 2, ksizet_l + 2, size5, size6/)
    subsizes = (/kthird, ksizex_l, ksizey_l, ksizet_l, size5, size6/)
    subsizes(direction + 2) = 1
    starts = (/0, 1, 1, 1, 0, 0/)
    starts(direction + 2) = position

    call MPI_Type_Create_Subarray(6, sizes, subsizes, starts, MPI_Order_Fortran, &
         & MPI_Double_Complex, typetarget, ierr)
    call MPI_Type_Commit(typetarget, ierr)
    return
  end subroutine init_single_halo_type_6

  subroutine init_halo_types_6()
    integer, parameter :: nsize5 = 1, nsize6 = 3
    integer :: size5(nsize5) = (/4/)
    integer :: size6(nsize6) = (/Nf, ndiag, ndiagg/)
    integer :: i5, i6

    do i5 = 1, nsize5
      do i6 = 1, nsize6
        call init_single_halo_type_6(0, 1, size5(i5), size6(i6), &
          & halo_6_xdn_send(size5(i5), size6(i6)))
        call init_single_halo_type_6(0, 0, size5(i5), size6(i6), &
          & halo_6_xdn_recv(size5(i5), size6(i6)))
        call init_single_halo_type_6(0, ksizex_l, size5(i5), size6(i6), &
          & halo_6_xup_send(size5(i5), size6(i6)))
        call init_single_halo_type_6(0, ksizex_l + 1, size5(i5), size6(i6), &
          & halo_6_xup_recv(size5(i5), size6(i6)))

        call init_single_halo_type_6(1, 1, size5(i5), size6(i6), &
          & halo_6_ydn_send(size5(i5), size6(i6)))
        call init_single_halo_type_6(1, 0, size5(i5), size6(i6), &
          & halo_6_ydn_recv(size5(i5), size6(i6)))
        call init_single_halo_type_6(1, ksizey_l, size5(i5), size6(i6), &
          & halo_6_yup_send(size5(i5), size6(i6)))
        call init_single_halo_type_6(1, ksizey_l + 1, size5(i5), size6(i6), &
          & halo_6_yup_recv(size5(i5), size6(i6)))

        call init_single_halo_type_6(2, 1, size5(i5), size6(i6), &
          & halo_6_tdn_send(size5(i5), size6(i6)))
        call init_single_halo_type_6(2, 0, size5(i5), size6(i6), &
          & halo_6_tdn_recv(size5(i5), size6(i6)))
        call init_single_halo_type_6(2, ksizet_l, size5(i5), size6(i6), &
          & halo_6_tup_send(size5(i5), size6(i6)))
        call init_single_halo_type_6(2, ksizet_l + 1, size5(i5), size6(i6), &
          & halo_6_tup_recv(size5(i5), size6(i6)))
      end do
    end do
  end subroutine init_halo_types_6
  !***********************************************************************
  !   Update boundary terms
  !***********************************************************************
  subroutine start_halo_update_6(size5, size6, Array, tag, reqs)
    !
    integer, intent(in) :: size5, size6, tag
    complex(dp), intent(inout) :: Array(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, &
        &                              0:ksizet_l + 1, size5, size6)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset
    integer :: ierr

    tag_offset = 6*tag

    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup, ierr)
    call MPI_Isend(Array, 1, halo_6_xup_send(size5, size6), ip_xup, 0 + tag_offset, comm, &
      & reqs(1), ierr)
    call MPI_Irecv(Array, 1, halo_6_xdn_recv(size5, size6), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2), ierr)
    call MPI_Isend(Array, 1, halo_6_xdn_send(size5, size6), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3), ierr)
    call MPI_Irecv(Array, 1, halo_6_xup_recv(size5, size6), ip_xup, 1 + tag_offset, comm, &
      & reqs(4), ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup, ierr)
    call MPI_Isend(Array, 1, halo_6_yup_send(size5, size6), ip_yup, 2 + tag_offset, comm, &
      & reqs(5), ierr)
    call MPI_Irecv(Array, 1, halo_6_ydn_recv(size5, size6), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6), ierr)
    call MPI_Isend(Array, 1, halo_6_ydn_send(size5, size6), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7), ierr)
    call MPI_Irecv(Array, 1, halo_6_yup_recv(size5, size6), ip_yup, 3 + tag_offset, comm, &
      & reqs(8), ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup, ierr)
    call MPI_Isend(Array, 1, halo_6_tup_send(size5, size6), ip_tup, 4 + tag_offset, comm, &
      & reqs(9), ierr)
    call MPI_Irecv(Array, 1, halo_6_tdn_recv(size5, size6), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10), ierr)
    call MPI_Isend(Array, 1, halo_6_tdn_send(size5, size6), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11), ierr)
    call MPI_Irecv(Array, 1, halo_6_tup_recv(size5, size6), ip_tup, 5 + tag_offset, comm, &
      & reqs(12), ierr)
  end subroutine start_halo_update_6

end module comms6
