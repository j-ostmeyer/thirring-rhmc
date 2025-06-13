module evolution
  implicit none

  real, parameter :: One = 1.0
contains

  subroutine initialise_phi_1flavour(phi, am, imass)
    use comms5, only: start_halo_update_5
    use counters, only: ancgpf, ancgpfpv
    use dum1, only: R, ps
    use mpi
    use gaussian, only: gauss0
    use params, only: kthird_l, ksizex_l, ksizey_l, ksizet_l, rescga, dp, ndiag
    use remez, only: anum4, aden4, bnum4, bden4
#ifdef HMC
    use measure_module, only: congrad
#else
    use qmrherm_module, only: qmrherm
#endif
    use dirac
    use vector
    use trial, only: u

    implicit none
    complex(dp), intent(out) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)!
    real, intent(in) :: am
    integer, intent(in) :: imass

    complex(dp) :: Xresult(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer :: reqs_ps(12)
    integer :: reqs_R(16)
    integer :: ierr
    integer :: idirac, ithird
    integer :: itercg
    !*******************************************************************
    !  Pseudofermion fields: Phi = {MdaggerM(1)}^-1/4 * {MdaggerM(m)}^1/4 * R, where
    !   R is gaussian
    !*******************************************************************

    do idirac = 1, 4
      do ithird = 1, kthird_l
#ifdef MPI
        call gauss0(ps, reqs_ps)
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)

#else
        call gauss0(ps)
#endif
        R(ithird, :, :, :, idirac) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))
      enddo
    enddo
    ! overkill because gauss0 should have already done the communications
    ! along the x,y and t axes, but we miss the "third" axis.
    call start_halo_update_5(4, R, 2342, reqs_R) ! "random" tag start
    !
#ifdef HMC
    call dslashd(Xresult,R,u,am,imass)
    call congrad(Xresult,rescga,itercg,One,1)
    ancgpfpv = ancgpfpv + float(itercg)
    call dslash(Phi,X,u,One,1)
#else
    !  For now Phi = {MdaggerM}^0.25 * R
    !
    call qmrherm(R, Xresult, rescga, itercg, am, imass, anum4, aden4, ndiag, 0)
    ancgpf = ancgpf + float(itercg)
    !
    R = Xresult
    !
    call qmrherm(R, Xresult, rescga, itercg, One, 1, bnum4, bden4, ndiag, 0)
    ancgpfpv = ancgpfpv + float(itercg)
    !
    Phi = Xresult
#endif

  end subroutine initialise_phi_1flavour

  subroutine initialise_pp(pp)
    use dum1, only: ps
    use comms, only: ip_third, comm_grp_third
    use gaussian, only: gaussp
    use mpi
    use params, only: ksizex_l, ksizey_l, ksizet_l
    implicit none
    real, intent(out) :: pp(ksizex_l, ksizey_l, ksizet_l, 3)
    integer :: reqs_ps(12)
    integer :: ierr
    integer :: mu
    !*******************************************************************
    !     heatbath for p
    !*******************************************************************
    !  for some occult reason this write statement is needed to ensure compatibility with earlier versions
    !     write(6,*) idum
    !     write(98,*) idum
    do mu = 1, 3
#ifdef MPI
      if (ip_third .eq. 0) then
        call gaussp(ps, reqs_ps)
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)
      endif
      ! we need to wait for the ranks with ip_third = 0 - Barrier is implicit
      ! We also need to broadcast halos, but we don't need the second component
      call MPI_Bcast(ps, &
                     (ksizex_l + 2)*(ksizey_l + 2)*(ksizet_l + 2), &
                     MPI_Real, 0, comm_grp_third, ierr)
#else
      call gaussp(ps)
#endif
      pp(:, :, :, mu) = ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1)
    enddo

  end subroutine initialise_pp

  subroutine evolve_theta_pp(iterl, dt, pp, phi, res, am, imass, itot)
    use force_module, only: force
    use gauge, only: coef
    use gforce, only: dSdpi
    use mpi
    use params, only: kthird_l, ksizex_l, ksizey_l, ksizet_l, dp
    use random
    ! we import these to be consistent with force, which imports them.
    use trial, only: ut => u, thetat => theta
    implicit none
    integer, intent(in) :: iterl
    real, intent(in) :: dt
    real, intent(inout) :: pp(ksizex_l, ksizey_l, ksizet_l, 3)
    complex(dp) :: phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(in) :: res
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(inout) :: itot

    integer :: iter, ierr
    real :: ytest, proby

    print *, rescgg

    !*******************************************************************
    !      half-step forward for p
    !*******************************************************************
    call force(Phi, res, am, imass)
    pp = pp - real(0.5*dt*dSdpi)
    !*******************************************************************
    !     start of main loop for classical time evolution
    !*******************************************************************
    proby = 1.0/float(iterl)

    do iter = 1, (4*iterl)
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
        write (6, "(A15,I4,A15,I4)") "MD iteration", iter, "of (average)", iterl
#ifdef MPI
      endif
#endif
      !  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
      thetat = thetat + dt*pp
      !  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
      call coef(ut, thetat)
      call force(phi, res, am, imass)

      ! test for end of random trajectory
      if (ip_global .eq. 0) then
        ytest = rano(yran, idum, 1, 1, 1)
      end if
#ifdef MPI
      call MPI_Bcast(ytest, 1, MPI_Real, 0, comm, ierr)
#endif
      if (ytest .lt. proby) then
        pp = pp - real(0.5*dt*dSdpi)
        exit
      else
        pp = pp - real(dt*dSdpi)
      endif

    enddo
    itot = itot + min(iter, 4*iterl)

  end subroutine evolve_theta_pp

end module evolution
