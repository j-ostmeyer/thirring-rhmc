module comms
  use params
  use pmpi
  implicit none

  save
  ! Fake arrays, just big enough to accomodate all pieces in a
  ! "start_halo_update" separately.
  integer, parameter :: fakesize4 = (ksizex_l+2)*(ksizey_l+2)*(ksizet_l+2)*4
  integer, parameter :: fakesize5 = kthird*fakesize4
  integer, parameter :: fakesize6 = kthird*fakesize4*ndiag
  real(dp) :: fakearray4r(fakesize4)
  complex(dp) :: fakearray4(fakesize4)
  complex(dp) :: fakearray5(fakesize5)
  complex(dp) :: fakearray6(fakesize6)

  integer :: ip_x, ip_y, ip_t, ip_global, np_global
  integer :: comm
  integer :: mpiio_type

  integer :: halo_6_xup_send(4:4, 1:25), halo_6_xdn_send(4:4, 1:25)
  integer :: halo_6_xup_recv(4:4, 1:25), halo_6_xdn_recv(4:4, 1:25)
  integer :: halo_6_yup_send(4:4, 1:25), halo_6_ydn_send(4:4, 1:25)
  integer :: halo_6_yup_recv(4:4, 1:25), halo_6_ydn_recv(4:4, 1:25)
  integer :: halo_6_tup_send(4:4, 1:25), halo_6_tdn_send(4:4, 1:25)
  integer :: halo_6_tup_recv(4:4, 1:25), halo_6_tdn_recv(4:4, 1:25)

  integer :: halo_5_xup_send(4:4), halo_5_xdn_send(4:4)
  integer :: halo_5_xup_recv(4:4), halo_5_xdn_recv(4:4)
  integer :: halo_5_yup_send(4:4), halo_5_ydn_send(4:4)
  integer :: halo_5_yup_recv(4:4), halo_5_ydn_recv(4:4)
  integer :: halo_5_tup_send(4:4), halo_5_tdn_send(4:4)
  integer :: halo_5_tup_recv(4:4), halo_5_tdn_recv(4:4)

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
  !***********************************************************************
  !   Define MPI datatypes to allow halo transfer
  !***********************************************************************
  subroutine init_halo_types()
    call init_halo_types_4
    call init_halo_types_4_real
    call init_halo_types_5
    call init_halo_types_6
  end subroutine init_halo_types

  ! Initialise a single halo type.
  ! direction: direction to be communicated; x=0, y=1, t=2
  ! position: coordinate in direction of the region
  ! i5, i6: size in extra directions
  ! datatype: datatype to transfer
  ! typetarget: array element to put type into;
  !   should be halo_{4,4_real,5,6}_{x,y,t}{up,dn}_{send,recv}(size{4,5}[, size6])
  subroutine init_single_halo_type_4(direction, position, size4, datatype, typetarget,offset)
    integer, intent(in) :: direction, position, size4
    integer, intent(in) :: datatype
    integer, intent(inout) :: offset
    integer, intent(out) :: typetarget
    integer :: ierr,idx
    integer :: size_to_transfer

    integer, dimension(4) :: subsizes
    subsizes = (/ ksizex_l, ksizey_l, ksizet_l, size4 /)
    subsizes(direction+1) = 1
    size_to_transfer = subsizes(1)
    do idx = 2,4
    size_to_transfer = size_to_transfer * subsizes(idx)
    enddo

    call MPI_Type_Create_Subarray(1, fakesize4, size_to_transfer, offset, MPI_Order_Fortran, &
      & datatype, typetarget,ierr)
    call MPI_Type_Commit(typetarget,ierr)
    offset = offset + size_to_transfer
    return
  end subroutine init_single_halo_type_4

  subroutine init_single_halo_type_5(direction, position, size5, typetarget,offset)
    integer, intent(in) :: direction, position, size5
    integer, intent(out) :: typetarget
    integer, intent(inout) :: offset
    integer :: ierr,idx
    integer :: size_to_transfer

    integer, dimension(5) :: subsizes
    subsizes = (/ kthird, ksizex_l, ksizey_l, ksizet_l, size5 /)
    subsizes(direction+1) = 1
    size_to_transfer = subsizes(1)
    do idx = 2,5
    size_to_transfer = size_to_transfer * subsizes(idx)
    enddo


    call MPI_Type_Create_Subarray(1, fakesize5, size_to_transfer, offset, MPI_Order_Fortran, &
      & MPI_Double_Complex, typetarget,ierr)
    call MPI_Type_Commit(typetarget,ierr)
    offset = offset + size_to_transfer
    return
  end subroutine init_single_halo_type_5

  subroutine init_single_halo_type_6(direction, position, size5, size6, typetarget,offset)
    integer, intent(in) :: direction, position, size5, size6
    integer, intent(out) :: typetarget
    integer, intent(inout) :: offset
    integer :: ierr,idx
    integer :: size_to_transfer

    integer, dimension(6) :: subsizes
    subsizes = (/ kthird, ksizex_l, ksizey_l, ksizet_l, size5, size6 /)
    subsizes(direction+1) = 1
    size_to_transfer = subsizes(1)
    do idx = 2,6
    size_to_transfer = size_to_transfer * subsizes(idx)
    enddo

    call MPI_Type_Create_Subarray(1, fakesize6, size_to_transfer, offset, MPI_Order_Fortran, &
      & MPI_Double_Complex, typetarget,ierr)
    call MPI_Type_Commit(typetarget,ierr)
    offset = offset + size_to_transfer
    return
  end subroutine init_single_halo_type_6

  subroutine init_halo_types_4()
    integer, parameter :: nsize4=2
    integer :: size4(nsize4) = (/ 3, 4 /)
    integer :: i4,offset
    integer :: type

    type = MPI_Double_Complex


    do i4 = 1,nsize4
      offset = 0
      call init_single_halo_type_4(0, 1,          size4(i4), type, &
        & halo_4_xdn_send(size4(i4)),offset)
      call init_single_halo_type_4(0, 0,          size4(i4), type, &
        & halo_4_xdn_recv(size4(i4)),offset)
      call init_single_halo_type_4(0, ksizex_l,   size4(i4), type, &
        & halo_4_xup_send(size4(i4)),offset)
      call init_single_halo_type_4(0, ksizex_l+1, size4(i4), type, &
        & halo_4_xup_recv(size4(i4)),offset)

      call init_single_halo_type_4(1, 1,          size4(i4), type, &
        & halo_4_ydn_send(size4(i4)),offset)
      call init_single_halo_type_4(1, 0,          size4(i4), type, &
        & halo_4_ydn_recv(size4(i4)),offset)
      call init_single_halo_type_4(1, ksizey_l,   size4(i4), type, &
        & halo_4_yup_send(size4(i4)),offset)
      call init_single_halo_type_4(1, ksizey_l+1, size4(i4), type, &
        & halo_4_yup_recv(size4(i4)),offset)

      call init_single_halo_type_4(2, 1,          size4(i4), type, &
        & halo_4_tdn_send(size4(i4)),offset)
      call init_single_halo_type_4(2, 0,          size4(i4), type, &
        & halo_4_tdn_recv(size4(i4)),offset)
      call init_single_halo_type_4(2, ksizet_l,   size4(i4), type, &
        & halo_4_tup_send(size4(i4)),offset)
      call init_single_halo_type_4(2, ksizet_l+1, size4(i4), type, &
        & halo_4_tup_recv(size4(i4)),offset)
    end do
  end subroutine init_halo_types_4

  subroutine init_halo_types_4_real()
    integer, parameter :: nsize4=1
    integer :: size4(nsize4) = (/ 2 /)
    integer :: i4,offset
    integer :: type

    type = MPI_Real

    do i4 = 1,nsize4
      offset = 0
      call init_single_halo_type_4(0, 1,          size4(i4), type, &
        & halo_4_real_xdn_send(size4(i4)),offset)
      call init_single_halo_type_4(0, 0,          size4(i4), type, &
        & halo_4_real_xdn_recv(size4(i4)),offset)
      call init_single_halo_type_4(0, ksizex_l,   size4(i4), type, &
        & halo_4_real_xup_send(size4(i4)),offset)
      call init_single_halo_type_4(0, ksizex_l+1, size4(i4), type, &
        & halo_4_real_xup_recv(size4(i4)),offset)

      call init_single_halo_type_4(1, 1,          size4(i4), type, &
        & halo_4_real_ydn_send(size4(i4)),offset)
      call init_single_halo_type_4(1, 0,          size4(i4), type, &
        & halo_4_real_ydn_recv(size4(i4)),offset)
      call init_single_halo_type_4(1, ksizey_l,   size4(i4), type, &
        & halo_4_real_yup_send(size4(i4)),offset)
      call init_single_halo_type_4(1, ksizey_l+1, size4(i4), type, &
        & halo_4_real_yup_recv(size4(i4)),offset)

      call init_single_halo_type_4(2, 1,          size4(i4), type, &
        & halo_4_real_tdn_send(size4(i4)),offset)
      call init_single_halo_type_4(2, 0,          size4(i4), type, &
        & halo_4_real_tdn_recv(size4(i4)),offset)
      call init_single_halo_type_4(2, ksizet_l,   size4(i4), type, &
        & halo_4_real_tup_send(size4(i4)),offset)
      call init_single_halo_type_4(2, ksizet_l+1, size4(i4), type, &
        & halo_4_real_tup_recv(size4(i4)),offset)
    end do
  end subroutine init_halo_types_4_real

  subroutine init_halo_types_5()
    integer, parameter :: nsize5=1
    integer :: size5(nsize5) = (/ 4 /)
    integer :: i5, offset

    do i5 = 1,nsize5
    offset = 0
    call init_single_halo_type_5(0, 1,          size5(i5), halo_5_xdn_send(size5(i5)),offset)
    call init_single_halo_type_5(0, 0,          size5(i5), halo_5_xdn_recv(size5(i5)),offset)
    call init_single_halo_type_5(0, ksizex_l,   size5(i5), halo_5_xup_send(size5(i5)),offset)
    call init_single_halo_type_5(0, ksizex_l+1, size5(i5), halo_5_xup_recv(size5(i5)),offset)

    call init_single_halo_type_5(1, 1,          size5(i5), halo_5_ydn_send(size5(i5)),offset)
    call init_single_halo_type_5(1, 0,          size5(i5), halo_5_ydn_recv(size5(i5)),offset)
    call init_single_halo_type_5(1, ksizey_l,   size5(i5), halo_5_yup_send(size5(i5)),offset)
    call init_single_halo_type_5(1, ksizey_l+1, size5(i5), halo_5_yup_recv(size5(i5)),offset)

    call init_single_halo_type_5(2, 1,          size5(i5), halo_5_tdn_send(size5(i5)),offset)
    call init_single_halo_type_5(2, 0,          size5(i5), halo_5_tdn_recv(size5(i5)),offset)
    call init_single_halo_type_5(2, ksizet_l,   size5(i5), halo_5_tup_send(size5(i5)),offset)
    call init_single_halo_type_5(2, ksizet_l+1, size5(i5), halo_5_tup_recv(size5(i5)),offset)
    end do
  end subroutine init_halo_types_5

  subroutine init_halo_types_6()
    integer, parameter :: nsize5=1, nsize6=3
    integer :: size5(nsize5) = (/ 4 /)
    integer :: size6(nsize6) = (/ Nf, ndiag, ndiagg /)
    integer :: i5, i6,offset

    do i5 = 1,nsize5
      do i6=1, nsize6
        offset = 0
        call init_single_halo_type_6(0, 1, size5(i5), size6(i6), &
          & halo_6_xdn_send(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(0, 0, size5(i5), size6(i6), &
          & halo_6_xdn_recv(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(0, ksizex_l, size5(i5), size6(i6), &
          & halo_6_xup_send(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(0, ksizex_l+1, size5(i5), size6(i6), &
          & halo_6_xup_recv(size5(i5), size6(i6)),offset)

        call init_single_halo_type_6(1, 1, size5(i5), size6(i6), &
          & halo_6_ydn_send(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(1, 0, size5(i5), size6(i6), &
          & halo_6_ydn_recv(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(1, ksizey_l, size5(i5), size6(i6), &
          & halo_6_yup_send(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(1, ksizey_l+1, size5(i5), size6(i6), &
          & halo_6_yup_recv(size5(i5), size6(i6)),offset)

        call init_single_halo_type_6(2, 1, size5(i5), size6(i6), &
          & halo_6_tdn_send(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(2, 0, size5(i5), size6(i6), &
          & halo_6_tdn_recv(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(2, ksizet_l, size5(i5), size6(i6), &
          & halo_6_tup_send(size5(i5), size6(i6)),offset)
        call init_single_halo_type_6(2, ksizet_l+1, size5(i5), size6(i6), &
          & halo_6_tup_recv(size5(i5), size6(i6)),offset)
      end do
    end do
  end subroutine init_halo_types_6
  !***********************************************************************
  !   Update boundary terms
  !***********************************************************************
  subroutine start_halo_update_4(size4, Array, tag, reqs)
    !     
    integer, intent(in) :: size4, tag
    complex(dp), intent(inout) :: Array(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, size4)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset 
    integer :: ierr

    tag_offset = 6 * tag

#ifdef SWAPSR
    !RECEIVES
    ! Start receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_xdn_recv(size4), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_xup_recv(size4), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_ydn_recv(size4), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_yup_recv(size4), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_tdn_recv(size4), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_tup_recv(size4), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)


    !SENDS
    ! Start send in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray4, 1, halo_4_xup_send(size4), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Isend(fakearray4, 1, halo_4_xdn_send(size4), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)

    ! Start send in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray4, 1, halo_4_yup_send(size4), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Isend(fakearray4, 1, halo_4_ydn_send(size4), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)

    ! Start send in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray4, 1, halo_4_tup_send(size4), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Isend(fakearray4, 1, halo_4_tdn_send(size4), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
    !      
#else
    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray4, 1, halo_4_xup_send(size4), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_xdn_recv(size4), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Isend(fakearray4, 1, halo_4_xdn_send(size4), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_xup_recv(size4), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray4, 1, halo_4_yup_send(size4), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_ydn_recv(size4), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Isend(fakearray4, 1, halo_4_ydn_send(size4), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_yup_recv(size4), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray4, 1, halo_4_tup_send(size4), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_tdn_recv(size4), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Isend(fakearray4, 1, halo_4_tdn_send(size4), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
    call MPI_Irecv(fakearray4, 1, halo_4_tup_recv(size4), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)

#endif
  end subroutine start_halo_update_4

  subroutine start_halo_update_4_real(size4, Array, tag, reqs)
    !     
    integer, intent(in) :: size4, tag
    real, intent(inout) :: Array(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, size4)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset 
    integer :: ierr

    tag_offset = 6 * tag

#ifdef SWAPSR
    ! RECEIVES
    ! Start receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_xdn_recv(size4), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_xup_recv(size4), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_ydn_recv(size4), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_yup_recv(size4), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_tdn_recv(size4), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_tup_recv(size4), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)

    ! SENDS
    ! Start send in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_xup_send(size4), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_xdn_send(size4), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)

    ! Start send in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_yup_send(size4), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_ydn_send(size4), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)

    ! Start send in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_tup_send(size4), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_tdn_send(size4), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
#else
    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_xup_send(size4), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_xdn_recv(size4), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_xdn_send(size4), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_xup_recv(size4), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_yup_send(size4), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_ydn_recv(size4), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_ydn_send(size4), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_yup_recv(size4), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_tup_send(size4), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_tdn_recv(size4), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Isend(fakearray4r, 1, halo_4_real_tdn_send(size4), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
    call MPI_Irecv(fakearray4r, 1, halo_4_real_tup_recv(size4), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)
#endif
    !      
  end subroutine start_halo_update_4_real

  subroutine start_halo_update_5(size5, Array, tag, reqs)
    !     
    integer, intent(in) :: size5, tag
    complex(dp), intent(inout) :: Array(kthird, 0:ksizex_l+1, 0:ksizey_l+1, &
      &                              0:ksizet_l+1, size5)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset 
    integer :: ierr

    tag_offset = 6 * tag

#ifdef SWAPSR
    ! RECEIVE
    ! Start receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_xdn_recv(size5), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_xup_recv(size5), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_ydn_recv(size5), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_yup_recv(size5), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_tdn_recv(size5), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_tup_recv(size5), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)

    !SENDS
    ! Start send in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray5, 1, halo_5_xup_send(size5), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Isend(fakearray5, 1, halo_5_xdn_send(size5), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)

    ! Start send in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray5, 1, halo_5_yup_send(size5), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Isend(fakearray5, 1, halo_5_ydn_send(size5), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)

    ! Start send in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray5, 1, halo_5_tup_send(size5), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Isend(fakearray5, 1, halo_5_tdn_send(size5), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
#else
    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray5, 1, halo_5_xup_send(size5), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_xdn_recv(size5), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Isend(fakearray5, 1, halo_5_xdn_send(size5), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_xup_recv(size5), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray5, 1, halo_5_yup_send(size5), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_ydn_recv(size5), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Isend(fakearray5, 1, halo_5_ydn_send(size5), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_yup_recv(size5), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray5, 1, halo_5_tup_send(size5), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_tdn_recv(size5), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Isend(fakearray5, 1, halo_5_tdn_send(size5), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
    call MPI_Irecv(fakearray5, 1, halo_5_tup_recv(size5), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)
#endif
    !      
  end subroutine start_halo_update_5

  subroutine start_halo_update_6(size5, size6, Array, tag, reqs)
    !     
    integer, intent(in) :: size5, size6, tag
    complex(dp), intent(inout) :: Array(kthird, 0:ksizex_l+1, 0:ksizey_l+1, &
      &                              0:ksizet_l+1, size5, size6)
    integer, intent(out) :: reqs(12)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn
    integer :: tag_offset 
    integer :: ierr

    tag_offset = 6 * tag

#ifdef SWAPSR
    ! RECEIVES
    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_xdn_recv(size5, size6), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_xup_recv(size5, size6), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_ydn_recv(size5, size6), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_yup_recv(size5, size6), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_tdn_recv(size5, size6), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_tup_recv(size5, size6), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)

    ! SENDS
    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray6, 1, halo_6_xup_send(size5, size6), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Isend(fakearray6, 1, halo_6_xdn_send(size5, size6), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray6, 1, halo_6_yup_send(size5, size6), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Isend(fakearray6, 1, halo_6_ydn_send(size5, size6), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray6, 1, halo_6_tup_send(size5, size6), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Isend(fakearray6, 1, halo_6_tdn_send(size5, size6), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
#else
    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup,ierr)
    call MPI_Isend(fakearray6, 1, halo_6_xup_send(size5, size6), ip_xup, 0 + tag_offset, comm, &
      & reqs(1),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_xdn_recv(size5, size6), ip_xdn, 0 + tag_offset, comm, &
      & reqs(2),ierr)
    call MPI_Isend(fakearray6, 1, halo_6_xdn_send(size5, size6), ip_xdn, 1 + tag_offset, comm, &
      & reqs(3),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_xup_recv(size5, size6), ip_xup, 1 + tag_offset, comm, &
      & reqs(4),ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup,ierr)
    call MPI_Isend(fakearray6, 1, halo_6_yup_send(size5, size6), ip_yup, 2 + tag_offset, comm, &
      & reqs(5),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_ydn_recv(size5, size6), ip_ydn, 2 + tag_offset, comm, &
      & reqs(6),ierr)
    call MPI_Isend(fakearray6, 1, halo_6_ydn_send(size5, size6), ip_ydn, 3 + tag_offset, comm, &
      & reqs(7),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_yup_recv(size5, size6), ip_yup, 3 + tag_offset, comm, &
      & reqs(8),ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup,ierr)
    call MPI_Isend(fakearray6, 1, halo_6_tup_send(size5, size6), ip_tup, 4 + tag_offset, comm, &
      & reqs(9),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_tdn_recv(size5, size6), ip_tdn, 4 + tag_offset, comm, &
      & reqs(10),ierr)
    call MPI_Isend(fakearray6, 1, halo_6_tdn_send(size5, size6), ip_tdn, 5 + tag_offset, comm, &
      & reqs(11),ierr)
    call MPI_Irecv(fakearray6, 1, halo_6_tup_recv(size5, size6), ip_tup, 5 + tag_offset, comm, &
      & reqs(12),ierr)
#endif
    !      
  end subroutine start_halo_update_6

  subroutine complete_halo_update(reqs)
    integer, intent(inout) :: reqs(12)
    integer :: ierr
    call MPI_Waitall(12, reqs, MPI_Statuses_Ignore,ierr)
  end subroutine complete_halo_update

  !***********************************************************************
  !   Initialise MPI variables
  !***********************************************************************
  subroutine init_MPI()
    integer :: coords(3)
#ifdef WITH_MUST
    integer, parameter :: must_rank = 1
#else
    integer, parameter :: must_rank = 0
#endif
    integer :: ierr

    call MPI_init(ierr)

    ! Check that we have the right number of processes
    call MPI_comm_size(MPI_COMM_WORLD, np_global,ierr)
    call MPI_comm_rank(MPI_COMM_WORLD, ip_global,ierr)
    if (np_global .ne. NP_X * NP_Y * NP_T + must_rank) then
      print *,"MPI dimensionality mismatch: ", NP_X, "*", NP_Y, "*", NP_T, "!=", np_global
      call MPI_finalize(ierr)
      call exit(2)
    end if

    ! Set up a Cartesian communicator; periodic boundaries, allow reordering
    call MPI_cart_create(MPI_COMM_WORLD, 3, (/ NP_X, NP_Y, NP_T /), &
      & (/ .true., .true., .true. /), .true., comm, ierr)

    ! Know where I am
    call MPI_cart_coords(comm, ip_global, 3, coords,ierr)
    ip_x = coords(1)
    ip_y = coords(2)
    ip_t = coords(3)

    ! Prepare file format for MPI-IO
    call MPI_Type_Create_Subarray(4, &! dimensionality
      (/ ksize, ksize, ksizet, 3 /), &! global volume
      (/ ksizex_l, ksizey_l, ksizet_l, 3 /), &! local volume
      (/ ip_x * ksizex_l, ip_y * ksizey_l, ip_t * ksizet_l, 0 /), &! start location
      MPI_Order_Fortran, &! array ordering
      MPI_Real, &! datatype to store
      mpiio_type,ierr ) ! type descriptor for this subarray type
    call MPI_Type_Commit(mpiio_type,ierr)

    ! Prepare all of the halo types
    call init_halo_types
    return
  end subroutine init_MPI
end module comms
