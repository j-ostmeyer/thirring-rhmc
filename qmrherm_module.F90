module qmrherm_module
  use params
  implicit none

  complex(dp) :: vtild(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: q(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: pm1(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: qm1(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: p(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: x1(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: x3(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

  complex(dp), save :: Phi0(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, ndiag)
  logical :: printall, we_have_warned

contains

  subroutine qmrherm(Phi, x, res, itercg, am, imass, anum, aden, ndiagq, iflag, use_sp, cg_returns)
    use comms_common, only: ip_global
    use comms, only: complete_halo_update, MPI_COMM_WORLD
#ifdef MPI
    use comms5, only: start_halo_update_5
    use comms6, only: start_halo_update_6
#else
    use comms5, only: update_halo_5
    use comms6, only: update_halo_6
#endif
    use derivs_module
    use dirac
    use gforce
    use multishift_module, only: multishift_solver, multishift_solver_sp
    use params
    use trial, only: u
    complex(dp), intent(in) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(out) :: x(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, intent(in) :: imass, ndiagq, iflag
    real(dp), intent(in) :: anum(0:ndiagq), aden(ndiagq)
    real, intent(in) :: res, am
    logical, intent(in), optional :: use_sp
    integer, intent(out) :: itercg
    integer, intent(out), optional :: cg_returns(ndiagq)
    logical :: force_dp = .false.
    integer :: cg_returns_tmp(ndiagq)
    real(dp) :: coeff
    integer :: idiag
#ifdef MPI
    integer, dimension(16) :: reqs_x2, reqs_Phi0, reqs_R, reqs_x
    integer :: ierr
#endif

    ! Set use_sp_flags_included (sp is not implemented for the wilson version)
#if defined(GENERATE_WITH_WILSON)
    force_dp = .true.
    if (.not. we_have_warned) then
      we_have_warned = .true.
      print *, "WARNING: Single precision is not supported in qmrherm when using the GENERATE_WITH_WILSON flag. Forcing &
      & double precision."
    endif
#endif


    if (ndiagq .gt. ndiag) then
      print *, 'The qmrherm_module module currently requires ndiagq be greater than ndiagg.'
      print *, 'Please adjust it and recompile.'
#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#endif
      call exit(1)
    endif

    if (mod(kthird_l, 4) .ne. 0) then
      print *, 'QMRHERM: kthird_l must be divisible by 4 (due to shift)'
#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#endif
      call exit(1)
    endif

    x = anum(0)*Phi
    if (present(use_sp) .and. use_sp .and. (.not. force_dp)) then
      call multishift_solver_sp(u, am, imass, ndiagq, aden, anum(1:ndiagq), x1, Phi, res, max_qmr_iters, itercg, cg_returns_tmp)
    else
      call multishift_solver(u, am, imass, ndiagq, aden, anum(1:ndiagq), x1, Phi, res, max_qmr_iters, itercg, cg_returns_tmp)
    endif

    if (present(cg_returns)) then
      cg_returns = cg_returns_tmp
    endif

    if (iflag .lt. 2) then
      !     Now evaluate solution x=(MdaggerM)^p * Phi
      do idiag = 1, ndiagq
        x(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
                x(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
              + anum(idiag)*x1(:, :, :, :, :, idiag)
      enddo
#ifdef MPI
      ! x is a saved module variable, so must be updated to avoid polluting the parent function
      ! could this in principle be moved outside the function so we don't do it unnecessarily?
      ! but in that case we wouldn't be able to hide the communications
      call start_halo_update_5(4, x, 3, reqs_x)
#else
      call update_halo_5(4, x)
#endif

      !  update phi0 block if required...
      if (iflag .eq. 1) then
        Phi0(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, 1:ndiagq) = &
                    x1(:, :, :, :, :, 1:ndiagq)
#ifdef MPI
        ! No way to hide communications here unfortunately
        ! In principle this could be better interleaved with the x update
        ! but that would add extra branching, and this section is messy enough already
        call start_halo_update_6(4, ndiagq, Phi0, 4, reqs_Phi0)
        call complete_halo_update(reqs_Phi0)
#else
        call update_halo_6(4, ndiagq, Phi0)
#endif
      endif! if(iflag.eq.1) then

#ifdef MPI
      call complete_halo_update(reqs_x)
#endif

    else ! if(iflag.lt.2)then

      do idiag = 1, ndiagq
        !  x2 = M*x1
        R(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(:, :, :, :, :, idiag)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, R, 5, reqs_R)
        call complete_halo_update(reqs_R)
#else
        call update_halo_5(4, R)
#endif

        ! Communication of x2 generated here can be hidden if iflag isn't 2, while R is updated
        call dslash(x2, R, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, x2, 6, reqs_x2)
#else
        call update_halo_5(4, x2)
#endif

        if (iflag .eq. 2) then
          coeff = anum(idiag)
#ifdef MPI
          call complete_halo_update(reqs_x2)
#endif
          call derivs(R, x2, coeff, 0, am, imass)

        else ! if(iflag.eq.2)then
          coeff = -anum(idiag)
          R = Phi0(:, :, :, :, :, idiag)
#ifdef MPI
          call complete_halo_update(reqs_x2)
#endif
          call derivs(R, x2, coeff, 0, am, imass)

          ! Communication of x2 generated here can be hidden while R is updated
          call dslash(x2, R, u, am, imass)
#ifdef MPI
          call start_halo_update_5(4, x2, 7, reqs_x2)
#else
          call update_halo_5(4, x2)
#endif

          R(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(:, :, :, :, :, idiag)
#ifdef MPI
          call start_halo_update_5(4, R, 8, reqs_R)
          call complete_halo_update(reqs_x2)
          call complete_halo_update(reqs_R)
#else
          call update_halo_5(4, R)
#endif
          call derivs(x2, R, coeff, 1, am, imass)
        endif ! if(iflag.eq.2)then
      enddo ! do idiag=1, ndiagq
    endif ! if(iflag.lt.2)then , else

    if (ip_global .eq. 0 .and. printall) then
      if (present(use_sp) .and. use_sp .and. (.not. force_dp)) then
        print *, "[SP] Qmrherm iterations,res:", itercg, res
      else
        print *, "[DP] Qmrherm iterations,res:", itercg, res
      endif
    endif
    return
  end subroutine qmrherm
  !**********************************************************************
  !  iflag = 0 : evaluates Rdagger*(Mdagger)'*x2
  !  iflag = 1 : evaluates Rdagger*(M)'*x2
  !**********************************************************************
end module qmrherm_module
