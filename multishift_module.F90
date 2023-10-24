#include "kernel.h"

#ifdef SCOREPINST
#include "scorep/SCOREP_User.inc"
#endif

module multishift_module
  implicit none

  integer, parameter :: print_every = 100000
contains

  ! From https://arxiv.org/abs/hep-lat/9612014
  ! Krylov space solvers for shifted linear systems, B. Jegerlehner, 1996
  subroutine multishift_solver(u, am, imass, ndiagq, aden, anum, output, input, res,&
    &maxcg, cg_return, cg_returns)
#if defined(NEWKERNEL) && defined(WILSONKERNEL)
    use diracWilson
#endif
#if defined(NEWKERNEL) && defined(SHAMIRKERNEL)
    use diracShamir
#endif
#ifndef NEWKERNEL
    use dirac
#endif
    use params
    use reductions
    use mpi
    use comms, only: complete_halo_update
    use comms_common, only: comm, ip_global
    use comms5, only: start_halo_update_5, init_halo_update_5
    ! subroutine parameters
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: ndiagq
    real(dp), intent(in) :: aden(ndiagq)
    real(dp), intent(in) :: anum(ndiagq)
    complex(dp) :: output(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
    complex(dp) :: input(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(in) :: res
    integer, intent(in) :: maxcg
    integer, intent(out) :: cg_return
    integer, intent(out), optional :: cg_returns(ndiagq)

    ! temporary variables - large vectors
    complex(dp) ::         r(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) ::         h(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) ::         s(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) ::         p(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: shiftferm(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)

    ! temporary variables - short vectors
    real(dp) :: zeta_i(ndiag)
    real(dp) :: zeta_ii(ndiag)
    real(dp) :: zeta_iii(ndiag)
    real(dp) :: omegas(ndiag)
    real(dp) :: gammas(ndiag)
    logical :: flags(ndiag) ! flag(i) = .true. : output(i) still has not converged.

    ! temporary variables - scalars
    real(dp) :: alpha, delta, lambda, omega, omega_save, gammag

    ! temporary variables - utilities
    integer :: ishift, maxishift, minishift
    real(dp) :: correction(ndiagq)
    real(dp) :: source_norm
#ifdef MPI
    integer, dimension(12) :: reqs_h, reqs_p
    integer :: ierr
    real(dp) :: dp_reduction ! DEBUG
#endif
    integer :: idirac, it, iy, ix, iz
    integer, parameter :: shift = 4

#ifdef SCOREPINST
    SCOREP_USER_REGION_DEFINE(dirac_op)
    SCOREP_USER_REGION_DEFINE(post)
#endif

    r = input ! vector
    p = r     ! vector
    delta = sum(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
#ifdef MPI
    call MPI_AllReduce(delta, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
    delta = dp_reduction
#endif

    source_norm = delta
    omega = 1.0d0
    alpha = 0.0d0
    output = 0.0d0 ! vector
    ! Setting up persistent communication requests
    call init_halo_update_5(4, h, 1, reqs_h)
    call init_halo_update_5(4, p, 2, reqs_p)

    correction = abs(anum/aden) ! ndiagq-long
    correction = correction/exp(sum(log(correction))/ndiagq) !

    do ishift = 1, ndiagq
      flags(ishift) = .true.
      shiftferm(:, :, :, :, :, ishift) = input(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
      zeta_i(ishift) = 1.0d0
      zeta_ii(ishift) = 1.0d0
      gammas(ishift) = 0.0d0
    enddo
    gammag = 0.0d0
    cg_return = 0
    maxishift = ndiagq
    minishift = 1

    do while ((maxishift .gt. (minishift - 1)) .and. (cg_return .lt. maxcg))
      cg_return = cg_return + 1

      ! DIRAC OPERATOR
      call dslash(h, p, u, am, imass)
#ifdef MPI
      call MPI_Startall(12, reqs_h, ierr)
      !call complete_halo_update(reqs_h) ! Now this call happens in dslashd
      call dslashd(s, h, u, am, imass, reqs_h)
#else
      call update_halo_5(4, h)
      call dslashd(s, h, u, am, imass)
#endif

      alpha = sum(real(conjg(p(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) &
                      &                *s(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))

#ifdef MPI
      call MPI_AllReduce(alpha, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
      alpha = dp_reduction
#endif

      omega_save = omega
      omega = -delta/alpha

      do ishift = minishift, maxishift
        zeta_iii(ishift) = (zeta_i(ishift)*zeta_ii(ishift)*omega_save)/ &
          & (omega*gammag*(zeta_i(ishift) - zeta_ii(ishift)) + &
          &   zeta_i(ishift)*omega_save*(1.0 - aden(ishift)*omega))
        omegas(ishift) = omega*zeta_iii(ishift)/zeta_ii(ishift)
      enddo

      r = r + omega*s
      lambda = sum(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
#ifdef MPI
      call MPI_AllReduce(lambda, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
      lambda = dp_reduction
#endif
      gammag = lambda/delta

      p = r + gammag*p
      call MPI_Startall(12, reqs_p, ierr)

      gammas(minishift:maxishift) = gammag*zeta_iii(minishift:maxishift)* &
                                    omegas(minishift:maxishift)/(zeta_ii(minishift:maxishift)*omega)

#ifdef SCOREPINST
      SCOREP_USER_REGION_BEGIN(post, 'post',&
        &SCOREP_USER_REGION_TYPE_COMMON)
#endif
      do ishift = minishift, maxishift
        do idirac = 1, 4
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                do iz = 1, kthird, shift
                  output(iz:iz + shift - 1, ix, iy, it, idirac, ishift) = &
                    &output(iz:iz + shift - 1, ix, iy, it, idirac, ishift) -&
                    &shiftferm(iz:iz + shift - 1, ix, iy, it, idirac, ishift)*&
                    &omegas(ishift)
                  shiftferm(iz:iz + shift - 1, ix, iy, it, idirac, ishift) = &
                    &shiftferm(iz:iz + shift - 1, ix, iy, it, idirac, ishift)*&
                   & gammas(ishift) +&
                   & r(iz:iz + shift - 1, ix, iy, it, idirac)*zeta_iii(ishift)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      if (present(cg_returns)) then
        cg_returns(minishift:maxishift) = cg_return
      endif

      !! Original code
      !do ishift=minishift,maxishift
      !  output(:,:,:,:,:,ishift) = output(:,:,:,:,:,ishift)-&
      !    &shiftferm(:,:,:,:,:,ishift)*omegas(ishift)
      !  shiftferm(:,:,:,:,:,ishift) = shiftferm(:,:,:,:,:,ishift) * gammas(ishift)+&
      !   & r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:) * zeta_iii(ishift)
      !enddo
#ifdef SCOREPINST
      SCOREP_USER_REGION_END(post)
#endif

!      test : block
!        use inverter_utils, only: dirac_op_shifted
!        ! CHECK : Real residual seems to decrease faster than the one computed...
!        complex(dp) :: xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!        complex(dp) :: xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!        complex(dp) :: check(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!#ifdef MPI
!        integer, dimension(12) :: reqs_xin
!#endif
!
!        real(dp) :: checksum
!
!        if(mod(cg_return,10)==0) then
!          do ishift =1,ndiagq
!            xin(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = output(:,:,:,:,:,ishift)
!#ifdef MPI
!            call start_halo_update_5(4, xin, 10, reqs_xin)
!            call complete_halo_update(reqs_xin)
!#else
!            call update_halo_5(4, xin)
!#endif
!            call dirac_op_shifted(xout,xin,u,am,imass,real(aden(ishift)))
!
!            check = i
!              & xout(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
!            checksum = sum(abs(check)**2)
!#ifdef MPI
!            call MPI_AllReduce(checksum, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
!            checksum = dp_reduction
!#endif
!            if(ip_global.eq.0) then
!              print*, ishift,sqrt(delta*zeta_i(ishift)**2),checksum,res,real(aden(ishift))
!            endif
!          enddo
!        endif
!      end block test

      maxishift = 0
      minishift = ndiagq + 1
      do ishift = 1, ndiagq
        if (flags(ishift)) then
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            if (mod(cg_return, print_every) .eq. 0) then
              write (6, '(E9.1E3)', advance="no") sqrt(delta*zeta_ii(ishift)**2)*correction(ishift)
            endif
#ifdef MPI
          endif
#endif
          if (sqrt(delta*zeta_ii(ishift)**2)*correction(ishift) < res) then
            flags(ishift) = .false.
          else
            maxishift = max(ishift, maxishift)
            minishift = min(ishift, minishift)
          endif
        else
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            if (mod(cg_return, print_every) .eq. 0) then
              write (6, '(A9)', advance="no") "-"
            endif
#ifdef MPI
          endif
#endif
        endif
      enddo
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
        if (mod(cg_return, print_every) .eq. 0) then
          write (6, *) ""
        endif
#ifdef MPI
      endif
#endif

      zeta_i(minishift:maxishift) = zeta_ii(minishift:maxishift)
      zeta_ii(minishift:maxishift) = zeta_iii(minishift:maxishift)

      delta = lambda

      call complete_halo_update(reqs_p)

      !if(ip_global.eq.0) then
      !   print*, "iteration", cg_return, "of", maxcg
      !endif
    enddo

  end subroutine multishift_solver

  ! From https://arxiv.org/abs/hep-lat/9612014
  ! Krylov space solvers for shifted linear systems, B. Jegerlehner, 1996
  ! Single precision version
#ifndef NEWKERNEL
  subroutine multishift_solver_sp(udp, am, imass, ndiagq, aden, anum, outputdp, inputdp, res,&
    &maxcg, cg_return, cg_returns)
    use dirac_sp, only: dslash_sp, dslashd_sp
    use params
    use comms_common, only: comm, ip_global
    use comms, only: complete_halo_update
    use mpi
    use comms5_sp, only: init_halo_update_5_sp
    use reductions, only: reduce_real_dp5d
    ! subroutine parameters
    complex(dp), intent(in) :: udp(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(sp) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: ndiagq
    real(dp), intent(in) :: aden(ndiagq)
    real(dp), intent(in) :: anum(ndiagq)
    complex(dp) :: outputdp(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
    complex(sp) :: output(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
    complex(dp) :: inputdp(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(sp) :: input(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(in) :: res
    integer, intent(in) :: maxcg
    integer, intent(out) :: cg_return
    integer, intent(out), optional :: cg_returns(ndiagq)

    ! temporary variables - large vectors
    complex(sp) ::         r(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(sp) ::         h(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(sp) ::         s(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(sp) ::         p(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(sp) :: shiftferm(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)

    ! temporary variables - short vectors
    real(dp) :: zeta_i(ndiag)
    real(dp) :: zeta_ii(ndiag)
    real(dp) :: zeta_iii(ndiag)
    real(dp) :: omegas(ndiag)
    real(dp) :: gammas(ndiag)
    logical :: flags(ndiag) ! flag(i) = .true. : output(i) still has not converged.

    ! temporary variables - scalars
    real(dp) :: alpha, delta, lambda, omega, omega_save, gammag

    ! temporary variables - utilities
    integer :: ishift, maxishift, minishift
    real(dp) :: correction(ndiagq)
    real(dp) :: source_norm
#ifdef MPI
    integer, dimension(12) :: reqs_h, reqs_p
    integer :: ierr
    real(dp) :: dp_reduction ! DEBUG
#endif
    integer :: idirac, it, iy, ix, iz
    integer, parameter :: shift = 4

#ifdef SCOREPINST
    SCOREP_USER_REGION_DEFINE(dirac_op)
    SCOREP_USER_REGION_DEFINE(post)
#endif

    print *,"sp multishift - stop"
    stop

    ! convert to single precision
    u = udp  ! type conversion
    input = inputdp ! type conversion

    r = input ! vector
    p = r     ! vector
    delta = reduce_real_dp5d(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
#ifdef MPI
    call MPI_AllReduce(delta, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
    delta = dp_reduction
#endif

    source_norm = delta
    omega = 1.0d0
    alpha = 0.0d0
    output = 0.0 ! vector
    ! Setting up persistent communication requests
    call init_halo_update_5_sp(4, h, 1, reqs_h)
    call init_halo_update_5_sp(4, p, 2, reqs_p)

    correction = abs(anum/aden) ! ndiagq-long
    correction = correction/exp(sum(log(correction))/ndiagq) !

    do ishift = 1, ndiagq
      flags(ishift) = .true.
      shiftferm(:, :, :, :, :, ishift) = input(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
      zeta_i(ishift) = 1.0d0
      zeta_ii(ishift) = 1.0d0
      gammas(ishift) = 0.0d0
    enddo
    gammag = 0.0d0
    cg_return = 0
    if (present(cg_returns)) then
      cg_returns = 0
    endif
    maxishift = ndiagq
    minishift = 1

    do while ((maxishift .gt. (minishift - 1)) .and. (cg_return .lt. maxcg))
      cg_return = cg_return + 1

      ! DIRAC OPERATOR
      call dslash_sp(h, p, u, am, imass)
#ifdef MPI
      call MPI_Startall(12, reqs_h, ierr)
      !call complete_halo_update(reqs_h) ! Now this call happens in dslashd
      call dslashd_sp(s, h, u, am, imass, reqs_h)
#else
      call update_halo_5_sp(4, h)
      call dslashd_sp(s, h, u, am, imass)
#endif

      alpha = reduce_real_dp5d(real(conjg(p(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) &
                                   &                *s(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))

#ifdef MPI
      call MPI_AllReduce(alpha, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
      alpha = dp_reduction
#endif

      omega_save = omega
      omega = -delta/alpha

      do ishift = minishift, maxishift
        zeta_iii(ishift) = (zeta_i(ishift)*zeta_ii(ishift)*omega_save)/ &
          & (omega*gammag*(zeta_i(ishift) - zeta_ii(ishift)) + &
          &   zeta_i(ishift)*omega_save*(1.0 - aden(ishift)*omega))
        omegas(ishift) = omega*zeta_iii(ishift)/zeta_ii(ishift)
      enddo

      r = r + real(omega)*s
      lambda = reduce_real_dp5d(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)
#ifdef MPI
      call MPI_AllReduce(lambda, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
      lambda = dp_reduction
#endif
      gammag = lambda/delta

      p = r + real(gammag)*p
      call MPI_Startall(12, reqs_p, ierr)

      gammas(minishift:maxishift) = gammag*zeta_iii(minishift:maxishift)* &
                                    omegas(minishift:maxishift)/(zeta_ii(minishift:maxishift)*omega)

#ifdef SCOREPINST
      SCOREP_USER_REGION_BEGIN(post, 'post',&
        &SCOREP_USER_REGION_TYPE_COMMON)
#endif
      do ishift = minishift, maxishift
        do idirac = 1, 4
          do it = 1, ksizet_l
            do iy = 1, ksizey_l
              do ix = 1, ksizex_l
                do iz = 1, kthird, shift
                  output(iz:iz + shift - 1, ix, iy, it, idirac, ishift) = &
                    &output(iz:iz + shift - 1, ix, iy, it, idirac, ishift) -&
                    &shiftferm(iz:iz + shift - 1, ix, iy, it, idirac, ishift)*&
                    &real(omegas(ishift))
                  shiftferm(iz:iz + shift - 1, ix, iy, it, idirac, ishift) = &
                    &shiftferm(iz:iz + shift - 1, ix, iy, it, idirac, ishift)*&
                   & real(gammas(ishift)) +&
              & r(iz:iz + shift - 1, ix, iy, it, idirac)*real(zeta_iii(ishift))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      if (present(cg_returns)) then
        cg_returns(minishift:maxishift) = cg_return
      endif

      !! Original code
      !do ishift=minishift,maxishift
      !  output(:,:,:,:,:,ishift) = output(:,:,:,:,:,ishift)-&
      !    &shiftferm(:,:,:,:,:,ishift)*omegas(ishift)
      !  shiftferm(:,:,:,:,:,ishift) = shiftferm(:,:,:,:,:,ishift) * gammas(ishift)+&
      !   & r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:) * zeta_iii(ishift)
      !enddo
#ifdef SCOREPINST
      SCOREP_USER_REGION_END(post)
#endif

      maxishift = 0
      minishift = ndiagq + 1
      do ishift = 1, ndiagq
        if (flags(ishift)) then
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            if (mod(cg_return, print_every) .eq. 0) then
              write (6, '(E9.1E3)', advance="no") sqrt(delta*zeta_ii(ishift)**2)*correction(ishift)
            endif
#ifdef MPI
          endif
#endif
          if (sqrt(delta*zeta_ii(ishift)**2)*correction(ishift) < res) then
            flags(ishift) = .false.
          else
            maxishift = max(ishift, maxishift)
            minishift = min(ishift, minishift)
          endif
        else
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            if (mod(cg_return, print_every) .eq. 0) then
              write (6, '(A9)', advance="no") "-"
            endif
#ifdef MPI
          endif
#endif
        endif
      enddo
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
        if (mod(cg_return, print_every) .eq. 0) then
          write (6, *) ""
        endif
#ifdef MPI
      endif
#endif

      zeta_i(minishift:maxishift) = zeta_ii(minishift:maxishift)
      zeta_ii(minishift:maxishift) = zeta_iii(minishift:maxishift)

      delta = lambda

      call complete_halo_update(reqs_p)

      !if(ip_global.eq.0) then
      !   print*, "iteration", cg_return, "of", maxcg
      !endif
    enddo

    outputdp = output

  end subroutine multishift_solver_sp
#endif

end module multishift_module
