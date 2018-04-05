program test_gauss0
      use dwf3d_lib
      use gauge
#ifdef MPI
      use mpi
#endif
      use comms
      use random
      implicit none

      real :: ps(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
      integer :: reqs(12), ix, iy, it, ix2, iy2, it2, i, j
      integer, dimension(4) :: duplicate_position1, duplicate_position2
      logical :: has_duplicates = .false.
      real :: y, sumps, maxps, minps

      seed = 1.0

! initialise MPI
#ifdef MPI
      call init_MPI
#endif
      
      call init_random(seed)
               
! call function
      call gauss0(ps, reqs)
#ifdef MPI
      call complete_halo_update(reqs)
#endif

! check output
      sumps = sum(ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
      maxps = maxval(ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
      minps = minval(ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
      outer: do it = 1, ksizet_l
         do it2 = 1, ksizet_l
            do iy = 1, ksizey_l
               do iy2 = 1, ksizey_l
                  do ix = 1, ksizex_l
                     do ix2 = 1, ksizex_l
                        do i = 1, 2
                           do j = 1, 2
                              if (ix .eq. ix2 .and. iy .eq. iy2 .and. it .eq. it2 &
                                   .and. i .eq. j) cycle
                              if (ps(ix, iy, it, i) .eq. ps(ix2, iy2, it2, j)) then
                                 has_duplicates = .true.
                                 duplicate_position1 = (/ ix, iy, it, i /)
                                 duplicate_position2 = (/ ix2, iy2, it2, j /)
                                 exit outer
                              endif
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do outer
#ifdef MPI
      call MPI_AllReduce(MPI_IN_PLACE, sumps, 1, MPI_REAL, MPI_SUM, comm, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, maxps, 1, MPI_REAL, MPI_MAX, comm, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, minps, 1, MPI_REAL, MPI_MIN, comm, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, has_duplicates, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
#endif
      if (ip_global .eq. 0) then
         if (sumps .lt. 14.460 .or. sumps .gt. 14.461) then
            print *, 'sum looks wrong:', sumps
         end if
         if (maxps .lt. 2.59536 .or. maxps .gt. 2.59537) then
            print *, 'max looks wrong:', maxps
         end if
         if (minps .lt. -2.441575 .or. minps .gt. -2.441565) then
            print *, 'min looks wrong:', minps
         end if
         if (has_duplicates) then
            print *, 'duplicate random numbers observed at:'
            write(*, '(4i3)') duplicate_position1
            write(*, '(4i3)') duplicate_position2
         end if
      end if
      ! print *, 'Seed: ', seed
#ifdef MPI
      call MPI_Finalize(ierr)
#endif
end program
