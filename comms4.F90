module comms4
  use comms_common
  use params
  use mpi
  implicit none

  integer :: halo_4_xup_send(3:4), halo_4_xdn_send(3:4)
  integer :: halo_4_xup_recv(3:4), halo_4_xdn_recv(3:4)
  integer :: halo_4_yup_send(3:4), halo_4_ydn_send(3:4)
  integer :: halo_4_yup_recv(3:4), halo_4_ydn_recv(3:4)
  integer :: halo_4_tup_send(3:4), halo_4_tdn_send(3:4)
  integer :: halo_4_tup_recv(3:4), halo_4_tdn_recv(3:4)

  integer :: halo_4_real_xup_send(2:2), halo_4_real_xdn_send(2:2)
  integer :: halo_4_real_xup_recv(2:2), halo_4_real_xdn_recv(2:2)
  integer :: halo_4_real_yup_send(2:2), halo_4_real_ydn_send(2:2)
  integer :: halo_4_real_yup_recv(2:2), halo_4_real_ydn_recv(2:2)
  integer :: halo_4_real_tup_send(2:2), halo_4_real_tdn_send(2:2)
  integer :: halo_4_real_tup_recv(2:2), halo_4_real_tdn_recv(2:2)

contains

  ! Initialise a single halo type.
  ! direction: direction to be communicated; x=0, y=1, t=2
  ! position: coordinate in direction of the region
  ! i5, i6: size in extra directions
  ! datatype: datatype to transfer
  ! typetarget: array element to put type into;
  !   should be halo_{4,4_real,5,6}_{x,y,t}{up,dn}_{send,recv}(size{4,5}[, size6])

  subroutine init_single_halo_type_4(direction, position, size4, datatype, typetarget)
    integer, intent(in) :: direction, position, size4
    integer, intent(in) :: datatype
    integer, intent(out) :: typetarget
    integer, dimension(4) :: sizes, subsizes, starts
    integer :: ierr

    sizes = (/ksizex_l + 2, ksizey_l + 2, ksizet_l + 2, size4/)
    subsizes = (/ksizex_l, ksizey_l, ksizet_l, size4/)
    subsizes(direction + 1) = 1
    starts = (/1, 1, 1, 0/)
    starts(direction + 1) = position

    call MPI_Type_Create_Subarray(4, sizes, subsizes, starts, MPI_Order_Fortran, &
         & datatype, typetarget, ierr)
    call MPI_Type_Commit(typetarget, ierr)
    return
  end subroutine init_single_halo_type_4

  subroutine init_halo_types_4()
    integer, parameter :: nsize4 = 2
    integer :: size4(nsize4) = (/3, 4/)
    integer :: i4
    integer :: type

    type = MPI_Double_Complex

    do i4 = 1, nsize4
      call init_single_halo_type_4(0, 1, size4(i4), type, &
        & halo_4_xdn_send(size4(i4)))
      call init_single_halo_type_4(0, 0, size4(i4), type, &
        & halo_4_xdn_recv(size4(i4)))
      call init_single_halo_type_4(0, ksizex_l, size4(i4), type, &
        & halo_4_xup_send(size4(i4)))
      call init_single_halo_type_4(0, ksizex_l + 1, size4(i4), type, &
        & halo_4_xup_recv(size4(i4)))

      call init_single_halo_type_4(1, 1, size4(i4), type, &
        & halo_4_ydn_send(size4(i4)))
      call init_single_halo_type_4(1, 0, size4(i4), type, &
        & halo_4_ydn_recv(size4(i4)))
      call init_single_halo_type_4(1, ksizey_l, size4(i4), type, &
        & halo_4_yup_send(size4(i4)))
      call init_single_halo_type_4(1, ksizey_l + 1, size4(i4), type, &
        & halo_4_yup_recv(size4(i4)))

      call init_single_halo_type_4(2, 1, size4(i4), type, &
        & halo_4_tdn_send(size4(i4)))
      call init_single_halo_type_4(2, 0, size4(i4), type, &
        & halo_4_tdn_recv(size4(i4)))
      call init_single_halo_type_4(2, ksizet_l, size4(i4), type, &
        & halo_4_tup_send(size4(i4)))
      call init_single_halo_type_4(2, ksizet_l + 1, size4(i4), type, &
        & halo_4_tup_recv(size4(i4)))
    end do
  end subroutine init_halo_types_4

  subroutine init_halo_types_4_real()
    integer, parameter :: nsize4 = 1
    integer :: size4(nsize4) = (/2/)
    integer :: i4
    integer :: type

    type = MPI_Real

    do i4 = 1, nsize4
      call init_single_halo_type_4(0, 1, size4(i4), type, &
        & halo_4_real_xdn_send(size4(i4)))
      call init_single_halo_type_4(0, 0, size4(i4), type, &
        & halo_4_real_xdn_recv(size4(i4)))
      call init_single_halo_type_4(0, ksizex_l, size4(i4), type, &
        & halo_4_real_xup_send(size4(i4)))
      call init_single_halo_type_4(0, ksizex_l + 1, size4(i4), type, &
        & halo_4_real_xup_recv(size4(i4)))

      call init_single_halo_type_4(1, 1, size4(i4), type, &
        & halo_4_real_ydn_send(size4(i4)))
      call init_single_halo_type_4(1, 0, size4(i4), type, &
        & halo_4_real_ydn_recv(size4(i4)))
      call init_single_halo_type_4(1, ksizey_l, size4(i4), type, &
        & halo_4_real_yup_send(size4(i4)))
      call init_single_halo_type_4(1, ksizey_l + 1, size4(i4), type, &
        & halo_4_real_yup_recv(size4(i4)))

      call init_single_halo_type_4(2, 1, size4(i4), type, &
        & halo_4_real_tdn_send(size4(i4)))
      call init_single_halo_type_4(2, 0, size4(i4), type, &
        & halo_4_real_tdn_recv(size4(i4)))
      call init_single_halo_type_4(2, ksizet_l, size4(i4), type, &
        & halo_4_real_tup_send(size4(i4)))
      call init_single_halo_type_4(2, ksizet_l + 1, size4(i4), type, &
        & halo_4_real_tup_recv(size4(i4)))
    end do
  end subroutine init_halo_types_4_real

  !***********************************************************************
  !   Update boundary terms
  !***********************************************************************

  subroutine start_halo_update_4(size4, Array, tag, reqs)
    !
    integer, intent(in) :: size4, tag
    complex(dp), intent(inout) :: Array(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, size4)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset
    integer :: ierr

    tag_offset = 6*tag

    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup, ierr)
    call MPI_Isend(Array, 1, halo_4_xup_send(size4), ip_xup, 0 + tag_offset, comm, &
      & reqs(1), ierr)
    call MPI_Irecv(Array, 1, halo_4_xdn_recv(size4), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2), ierr)
    call MPI_Isend(Array, 1, halo_4_xdn_send(size4), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3), ierr)
    call MPI_Irecv(Array, 1, halo_4_xup_recv(size4), ip_xup, 1 + tag_offset, comm, &
      & reqs(4), ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup, ierr)
    call MPI_Isend(Array, 1, halo_4_yup_send(size4), ip_yup, 2 + tag_offset, comm, &
      & reqs(5), ierr)
    call MPI_Irecv(Array, 1, halo_4_ydn_recv(size4), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6), ierr)
    call MPI_Isend(Array, 1, halo_4_ydn_send(size4), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7), ierr)
    call MPI_Irecv(Array, 1, halo_4_yup_recv(size4), ip_yup, 3 + tag_offset, comm, &
      & reqs(8), ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup, ierr)
    call MPI_Isend(Array, 1, halo_4_tup_send(size4), ip_tup, 4 + tag_offset, comm, &
      & reqs(9), ierr)
    call MPI_Irecv(Array, 1, halo_4_tdn_recv(size4), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10), ierr)
    call MPI_Isend(Array, 1, halo_4_tdn_send(size4), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11), ierr)
    call MPI_Irecv(Array, 1, halo_4_tup_recv(size4), ip_tup, 5 + tag_offset, comm, &
      & reqs(12), ierr)

  end subroutine start_halo_update_4

  subroutine start_halo_update_4_real(size4, Array, tag, reqs)
    !
    integer, intent(in) :: size4, tag
    real, intent(inout) :: Array(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, size4)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset
    integer :: ierr

    tag_offset = 6*tag

    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup, ierr)
    call MPI_Isend(Array, 1, halo_4_real_xup_send(size4), ip_xup, 0 + tag_offset, comm, &
      & reqs(1), ierr)
    call MPI_Irecv(Array, 1, halo_4_real_xdn_recv(size4), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2), ierr)
    call MPI_Isend(Array, 1, halo_4_real_xdn_send(size4), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3), ierr)
    call MPI_Irecv(Array, 1, halo_4_real_xup_recv(size4), ip_xup, 1 + tag_offset, comm, &
      & reqs(4), ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup, ierr)
    call MPI_Isend(Array, 1, halo_4_real_yup_send(size4), ip_yup, 2 + tag_offset, comm, &
      & reqs(5), ierr)
    call MPI_Irecv(Array, 1, halo_4_real_ydn_recv(size4), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6), ierr)
    call MPI_Isend(Array, 1, halo_4_real_ydn_send(size4), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7), ierr)
    call MPI_Irecv(Array, 1, halo_4_real_yup_recv(size4), ip_yup, 3 + tag_offset, comm, &
      & reqs(8), ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup, ierr)
    call MPI_Isend(Array, 1, halo_4_real_tup_send(size4), ip_tup, 4 + tag_offset, comm, &
      & reqs(9), ierr)
    call MPI_Irecv(Array, 1, halo_4_real_tdn_recv(size4), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10), ierr)
    call MPI_Isend(Array, 1, halo_4_real_tdn_send(size4), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11), ierr)
    call MPI_Irecv(Array, 1, halo_4_real_tup_recv(size4), ip_tup, 5 + tag_offset, comm, &
      & reqs(12), ierr)
  end subroutine start_halo_update_4_real

end module comms4
