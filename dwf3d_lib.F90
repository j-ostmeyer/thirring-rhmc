module dwf3d_lib
  use params
  implicit none
  save

  ! Random numbers
  real(dp) :: seed

! Useful constants
  real, parameter :: One = 1.0
contains

  subroutine dwf3d_main
    use random
    use gaussian
    use remez
    use remezg
    use trial, ut => u, thetat => theta
    use gauge
    use vector, X1 => X
    use gforce
    use avgitercounts
    use dum1
    use comms
    use measure_module
    use qmrherm_module, only: qmrherm, qmrhprint => printall
!*******************************************************************
!    Rational Hybrid Monte Carlo algorithm for bulk Thirring Model with Domain Wall
!         fermions
!
!    requires operation of QMR on complex vectors to determine
!    (Mdagger M)**-1  Phi
!
!    requires input partial fraction coefficients from Remez algorithm
!    of Clarke & Kennedy
!
!    the "third" direction is actually indexed 4 in the code - sorry!
!      (should have used C I know)
!
!           { 1 - hermition mass term psibar psi
!    imass= { 3 - antih.    mass term i psibar gamma_3 psi
!           { 5 - antih.    mass term i psibar gamma_5 psi
!     linear combinations of above require code modification
!
!    code exploits fact that gamma_3 is diagonal in Dirac basis used - speeds up
!    evolution slightly, and measurement by factor of two.
!
!    Fermion expectation values are measured using a noisy estimator.
!    on the Wilson matrix, which has dimension 4*kvol*kthird
!    inversions done using congrad, and matrix multiplies with dslash, dslashd
!
!    Pauli-Villars determinant defined using a hermitian mass m_h=One
!
!    trajectory length is random with mean dt*iterl
!    The code runs for a fixed number iter2 of trajectories.
!
!    Phi: pseudofermion field
!    am: bare fermion mass
!    actiona: running average of total action
!
!                                               SJH February 2017
!*******************************************************************
    implicit none
!     complex qq,qbqb
!     complex u
!     complex a,b
    complex(dp) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)!
    complex(dp) :: Xresult(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
! Currently unused
!    complex(dp) :: qq,qbqb
!    complex(dp) :: u
!    complex(dp) :: a,b
    real(dp) :: H0, H1, S0, S1, dH, dS, hg, hp
    real :: action, paction, gaction
    real :: vel2, x, ytest, atraj
    real :: dt, am, y, traj, proby
    real :: actiona, vel2a, pbp, pbpa, yav, yyav
    real :: ancgm, ancgma
    integer :: imass, iter, iterl, iter2, i, ia, idirac, ithird
    integer :: naccp, ipbp, itot, isweep, itercg, mu
    integer ::walltimesec
!*******************************************************************
!     variables to keep track of MPI requests
!*******************************************************************
#ifdef MPI
    integer :: reqs_ps(12)
    integer :: ierr
#endif
!
!*******************************************************************
!     input
!*******************************************************************
    complex(dp), parameter :: zi = (0.0, 1.0)
    ibound = -1
    qmrhprint = .true.
#ifdef MPI
    call init_MPI
#endif
!*******************************************************************
!     end of input
!*******************************************************************
!*******************************************************************
!     check qmrherm is going to be OK
!*******************************************************************
    if (ndiagg .gt. ndiag) then
      print *, 'The qmrherm_module module currently requires ndiag be greater than ndiagg.'
      print *, 'Please adjust it and recompile.'
      call exit(1)
    endif
    if (ip_global .eq. 0) then
      open (unit=7, file='output', status='unknown')
      open (unit=98, file='control', status='unknown')
    end if
    open (unit=25, file='midout', status='old')
    open (unit=36, file='remez2', status='old')
    open (unit=37, file='remez4', status='old')
    open (unit=38, file='remez2g', status='old')
    open (unit=39, file='remez4g', status='old')
    if (iread .eq. 1) then
      call sread
    endif
    read (25, *) dt, beta, am3, am, imass, iterl, iter2, walltimesec
    close (25)
! set a new seed by hand...

    if (iseed .ne. 0) then
      seed = 4139764973254.0 ! this should be a double precision number
      ! ending with d0, but we do not care
    endif
    call readseed(seed)
    if (ip_global .eq. 0) then
      write (7, *) 'seed: ', seed
    end if
    call init_random(seed)

    if (ip_global .eq. 0) then
      write (6, *) 'Initialized, ', seed
    endif
!*******************************************************************
!     initialization
!     istart.lt.0 : start from tape
!     istart=0    : ordered start
!     istart=1    : random start
!*******************************************************************
    call init(istart)
!  read in Remez coefficients
    read (36, *) anum2(0)
    read (37, *) anum4(0)
    read (38, *) anum2g(0)
    read (39, *) anum4g(0)
    do i = 1, ndiag
      read (36, *) anum2(i), aden2(i)
      read (37, *) anum4(i), aden4(i)
    enddo
    do i = 1, ndiagg
      read (38, *) anum2g(i), aden2g(i)
      read (39, *) anum4g(i), aden4g(i)
    enddo
    read (36, *) bnum2(0)
    read (37, *) bnum4(0)
    read (38, *) bnum2g(0)
    read (39, *) bnum4g(0)
    do i = 1, ndiag
      read (36, *) bnum2(i), bden2(i)
      read (37, *) bnum4(i), bden4(i)
    enddo
    do i = 1, ndiagg
      read (38, *) bnum2g(i), bden2g(i)
      read (39, *) bnum4g(i), bden4g(i)
    enddo
!*******************************************************************
!     print heading
!*******************************************************************
    traj = iterl*dt
    proby = 1.0/float(iterl)
    if (ip_global .eq. 0) then
      write (7, 9001) ksize, ksizet, kthird, Nf, dt, traj, ndiag, ndiagg, &
           & iter2, beta, am3, am, imass
9001  format(' ksize=', i3, ' ksizet=', i3, / &
             , ' kthird=', i3, / &
             , ' Nf =', i3, / &
             , ' time step: dt=', f6.4, ' trajectory length=', f9.6, / &
             , ' Remez ndiag: action =', i3, ' guidance=', i3, / &
             , ' # trajectories=', i6, ' beta=', f9.6, / &
             , ' am3=', f6.4, ' am=', f6.4/ &
             , ' imass=', i2)
#ifdef MPI
      write (7, 9002) NP_X, NP_Y, NP_T, ksizex_l, ksizey_l, ksizet_l
9002  format(" NP_X=", i3, " NP_Y=", i3, " NP_T=", i3, / &
           & " ksizex_l=", i3, " ksizey_l=", i3, " ksizet_l=", i3)
#endif
!     write(6,9004) rescgg,rescga,respbp
      write (7, 9004) rescgg, rescga, respbp
9004  format(' Stopping residuals: guidance: ', e11.4, ' acceptance: ', &
           &     e11.4, ' estimator: ', e11.4)
!     write(6,9044) rescgm
      write (7, 9044) rescgm
9044  format(' Stopping residuals: meson: ', e11.4)
      call rranget(seed, 1, 1, 1)
! c     write(6,*) 'seed: ', seed
      write (7, *) 'seed: ', seed
    end if
!*******************************************************************
!       initialize for averages
!*******************************************************************
    actiona = 0.0
    vel2a = 0.0
    pbpa = 0.0
    ancg = 0.0
    ancgh = 0.0
    ancgf = 0.0
    ancgpf = 0.0
    ancgpv = 0.0
    ancghpv = 0.0
    ancgfpv = 0.0
    ancgpfpv = 0.0
    ancgma = 0.0
    yav = 0.0
    yyav = 0.0
    naccp = 0
    ipbp = 0
    itot = 0
!*******************************************************************
!     start of classical evolution
!*******************************************************************
    classical_evolution: block
      real(dp) :: run_time, time_per_md_step ! conservative estimates
      real(dp) :: measurement_time, total_md_time
      total_md_time = 0

      do isweep = 1, iter2

#ifdef MPI
        if (ip_global .eq. 0) then
#endif
          write (6, *) 'Isweep', isweep, ' of', iter2
#ifdef MPI
        endif
#endif
! uncomment line below to go straight to measurement
!     goto 666
!*******************************************************************
!     initialise trial fields
!*******************************************************************
        thetat = theta
!
        call coef(ut, thetat)
!*******************************************************************
!  Pseudofermion fields: Phi = {MdaggerM(1)}^-1/4 * {MdaggerM(m)}^1/4 * R, where
!   R is gaussian
!*******************************************************************
        do ia = 1, Nf
!
          do idirac = 1, 4
            do ithird = 1, kthird
#ifdef MPI
              call gauss0(ps, reqs_ps)
              call complete_halo_update(reqs_ps)
#else
              call gauss0(ps)
#endif
              R(ithird, :, :, :, idirac) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))
            enddo
          enddo
!
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
!
        enddo
!*******************************************************************
!     heatbath for p
!*******************************************************************
!  for some occult reason this write statement is needed to ensure compatibility with earlier versions
!     write(6,*) idum
!     write(98,*) idum
        do mu = 1, 3
#ifdef MPI
          call gaussp(ps, reqs_ps)
          call complete_halo_update(reqs_ps)
#else
          call gaussp(ps)
#endif
          pp(:, :, :, mu) = ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1)
        enddo

        ! Start computing time (including hamiltonian calls)
        call cpu_time(run_time)
        total_md_time = total_md_time - run_time

!*******************************************************************
!  call to Hamiltonian
!
        call hamilton(Phi, H0, hg, hp, S0, rescga, isweep, 0, am, imass)
        if (isweep .eq. 1) then
          action = real(S0)/kvol
          gaction = real(hg)/kvol
          paction = real(hp)/kvol
        endif
!     goto 501
!*******************************************************************
!      half-step forward for p
!*******************************************************************
        call force(Phi, rescgg, am, imass, isweep, 0)
        pp = pp - 0.5*dt*dSdpi
!*******************************************************************
!     start of main loop for classical time evolution
!*******************************************************************

        do iter = 1, (4*iterl)
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            write (6, "(A15,I4,A15,I4)") "MD iteration", iter, "of (average)", iterl
#ifdef MPI
          endif
#endif
!
!  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
!
          thetat = thetat + dt*pp
!
!  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
!
          call coef(ut, thetat)
          call force(Phi, rescgg, am, imass, isweep, iter)
!
! test for end of random trajectory
!
          if (ip_global .eq. 0) then
            ytest = rano(yran, idum, 1, 1, 1)
          end if
#ifdef MPI
          call MPI_Bcast(ytest, 1, MPI_Real, 0, comm, ierr)
#endif
          if (ytest .lt. proby) then
            pp = pp - 0.5*dt*dSdpi
            itot = itot + iter
            goto 501
          else
            pp = pp - dt*dSdpi
          endif
!
        end do
!**********************************************************************
!  Monte Carlo step: accept new fields with probability=
!              min(1,exp(H0-H1))
!**********************************************************************
501     continue

        call hamilton(Phi, H1, hg, hp, S1, rescga, isweep, -1, am, imass)
        dH = H0 - H1
        dS = S0 - S1
        if (ip_global .eq. 0) then
          write (98, *) dH, dS
          write (6, *) "dH,dS ", dH, dS
        end if
        y = exp(real(dH))
        yav = yav + y
        yyav = yyav + y*y
!
        if (dH .lt. 0.0) then
          x = rano(yran, idum, 1, 1, 1)
#ifdef MPI
          call MPI_Bcast(x, 1, MPI_Real, 0, comm, ierr)
#endif
          if (x .gt. y) goto 600
        endif
!
!     step accepted: set s=st
!
        theta = thetat
        naccp = naccp + 1
        action = real(S1)/kvol
        gaction = real(hg)/kvol
        paction = real(hp)/kvol
600     continue
        if (ip_global .eq. 0) then
          write (11, *) isweep, gaction, paction
        end if
        actiona = actiona + action
        vel2 = sum(pp*pp)
#ifdef MPI
        call MPI_AllReduce(MPI_In_Place, vel2, 1, MPI_Real, MPI_Sum, comm, ierr)
#endif
        vel2 = vel2/(3*kvol)
        vel2a = vel2a + vel2

        ! Including also hamiltonian call time
        call cpu_time(run_time)
        total_md_time = total_md_time + run_time
        time_per_md_step = total_md_time/itot

!     uncomment to disable measurements
!     goto 601
!666    continue

        if ((isweep/iprint)*iprint .eq. isweep) then
          call cpu_time(run_time)
          thetat = theta
          call coef(ut, thetat)
          call measure(pbp, respbp, ancgm, am, imass)
!         call meson(rescgm,itercg,ancgm,am,imass)
          pbpa = pbpa + pbp
          ancgma = ancgma + ancgm
          ipbp = ipbp + 1
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            write (6, *) isweep, 'pbp:', pbp, ancgm
#ifdef MPI
          endif
#endif
          call cpu_time(measurement_time)
          measurement_time = measurement_time - run_time
        endif
!
        if ((isweep/icheck)*icheck .eq. isweep) then
          call rranget(seed, 1, 1, 1)
          if (iwrite .eq. 1) then
            call swrite
          endif
          flush (100)
          flush (200)
        endif

#ifdef MPI
        call MPI_AllReduce(MPI_In_Place, measurement_time, 1, MPI_Double_Precision, MPI_Max, comm, ierr)
        call MPI_AllReduce(MPI_In_Place, time_per_md_step, 1, MPI_Double_Precision, MPI_Max, comm, ierr)
#endif
        keep_running_check: block
          real(dp) :: run_time_left, time_for_next_iteration

          time_for_next_iteration = time_per_md_step*4*iterl*2

          if (((isweep + 1)/icheck)*icheck .eq. isweep) then
            time_for_next_iteration = time_for_next_iteration + measurement_time
          endif

          call cpu_time(run_time)

          if (run_time + time_for_next_iteration .gt. walltimesec) then
            if(ip_global.eq.0)then
              print*,'Expected next run time:',&
               run_time + time_for_next_iteration, ' larger than ',&
               walltimesec
              print*,'Quitting'
            endif
            exit
          endif
        end block keep_running_check
      end do
    end block classical_evolution
!*******************************************************************
!     end of main loop
!*******************************************************************
    actiona = actiona/iter2
    vel2a = vel2a/iter2
    pbpa = pbpa/ipbp
    ancg = ancg/(Nf*itot)
    ancgh = ancgh/(2*Nf*iter2)
    ancgpf = ancgpf/(Nf*iter2)
    ancgpv = ancgpv/(Nf*itot)
    ancgf = ancgf/(Nf*itot)
    ancgfpv = ancgfpv/(Nf*itot)
    ancghpv = ancghpv/(2*Nf*iter2)
    ancgpfpv = ancgpfpv/(iter2*Nf)
    ancgma = ancgma/ipbp
    yav = yav/iter2
    yyav = yyav/iter2 - yav*yav
    yyav = sqrt(yyav/(iter2 - 1))
    atraj = dt*itot/iter2
!*******************************************************************
!     print global averages
!*******************************************************************
    if (ip_global .eq. 0) then
!     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgpv,ancgh,ancghpv,ancgf,
!    & ancgfpv,ancgpf,ancgpfpv,pbpa,vel2a,actiona
      write (7, 9022) iter2, naccp, atraj, yav, yyav, &
           & ancg, ancgpv, ancgh, ancghpv, ancgf, ancgfpv, ancgpf, ancgpfpv, &
           & pbpa, ancgma, vel2a, actiona
9022  format(' averages for last ', i6, ' trajectories', /  &
       & 1x, ' # of acceptances: ', i6, ' average trajectory length= ', f8.3/ &
           & 1x, ' <exp-dH>=', e11.4, ' +/-', e10.3/ &
           & 1x, ' av. # QMR itr.'/ &
           & 1x, '     guidance: DWF  ', f9.3, '; PV  ', f9.3/ &
           & 1x, '   acceptance: DWF  ', f9.3, '; PV  ', f9.3/ &
           & 1x, '        force: DWF  ', f9.3, '; PV  ', f9.3/ &
           & 1x, 'pseudofermion: DWF  ', f9.3, '; PV  ', f9.3/ &
           & 1x, ' psibarpsi=', e11.3/ &
           & 1x, ' av. # QMR itr.', f9.3// &
          & 1x, ' mean square velocity=', e10.3, '; action per site=', e10.3//)
      write (7, 9024)
      write (7, 9024)
9024  format(1x)
!
      close (11)
    end if
!
    if (iwrite .eq. 1) then
      call rranget(seed, 1, 1, 1)
      call swrite
      call saveseed(seed)
      write (7, *) 'seed: ', idum
    endif

#ifdef MPI
    call MPI_Finalize(ierr)
#endif
!
  end subroutine dwf3d_main
!******************************************************************
!   calculate dSds for gauge fields at each intermediate time
!******************************************************************
  subroutine force(Phi, res1, am, imass, isweep, iter)
    use remezg
    use trial
    use gforce
    use avgitercounts
    use comms
    use qmrherm_module, only: qmrherm

    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, Nf)
    real, intent(in) :: res1, am
    integer, intent(in) :: imass, isweep, iter
!     complex Phi(kferm,Nf),X2(kferm)
!     complex X1,u
    complex(dp) :: X2(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Xresult(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer :: ia, itercg
!
!     write(6,111)
!111 format(' Hi from force')
!
    dSdpi = 0.0
!
! uncomment this line to quench the fermions!
!     return
! pseudofermion action is
!   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
!
    do ia = 1, Nf
!
      X2 = Phi(:, :, :, :, :, ia)

      call qmrherm(X2, Xresult, res1, itercg, One, 1, anum4g, aden4g, ndiagg, 1, spmd)
      ancgpv = ancgpv + float(itercg)

      X2 = Xresult
!
      call qmrherm(X2, Xresult, res1, itercg, am, imass, bnum2g, bden2g, ndiagg, 0, spmd)
      ancg = ancg + float(itercg)
!     write(111,*) itercg
      X2 = Xresult
!
!  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
      call qmrherm(X2, Xresult, res1, itercg, am, imass, anum2g, aden2g, ndiagg, 2, spmd)
      ancgf = ancgf + float(itercg)

!     write(113,*) itercg
!  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
      call qmrherm(X2, Xresult, res1, itercg, One, 1, anum4g, aden4g, ndiagg, 3, spmd)
      ancgfpv = ancgfpv + float(itercg)
!
    enddo
!
    if (ibound .eq. -1 .and. ip_t .eq. (np_t - 1)) then
      dSdpi(:, :, ksizet_l, 3) = -dSdpi(:, :, ksizet_l, 3)
    endif

    dSdpi = dSdpi + beta*Nf*theta
!
    return
  end subroutine force
!******************************************************************
!   Evaluation of Hamiltonian function
!******************************************************************
  subroutine hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)
    use remez
    use trial, only: theta, pp
    use dum1
    use avgitercounts
    use comms
    use qmrherm_module, only: qmrherm
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, Nf)
    real(dp), intent(out) :: h, hg, hp, s
    real, intent(in) :: res2, am
    integer, intent(in) :: isweep, iflag, imass
    complex(dp) :: Xresult(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp) :: hf
    integer :: itercg, ia
#ifdef MPI
    integer :: ierr
#endif
!     write(6,111)
!111 format(' Hi from hamilton')
!
    hf = 0.0
!
    hp = 0.5*sum(pp**2)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, hp, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
#endif

    hg = 0.5*Nf*beta*sum(theta**2)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, hg, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
#endif
    h = hg + hp

!
! uncomment these lines to quench the fermions!
!     return
!
!  pseudofermion action is
!   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
!
    do ia = 1, Nf
!
      R = Phi(:, :, :, :, :, ia)

      call qmrherm(R, Xresult, res2, itercg, One, 1, anum4, aden4, ndiag, 0)
      ancghpv = ancghpv + float(itercg)
!
      R = Xresult
!
      call qmrherm(R, Xresult, res2, itercg, am, imass, bnum2, bden2, ndiag, 0)
      ancgh = ancgh + float(itercg)
!
      hf = hf + sum(real(conjg(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) &
                        &        *Xresult(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))
!
    enddo
#ifdef MPI
! hf is built up from zero during the loop so only needs to be summed across
! all partitions at this point
    call MPI_AllReduce(MPI_In_Place, hf, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
#endif
!
    h = hg + hp + hf
    s = hg + hf
!
    return
  end subroutine hamilton

! c******************************************************************
! c   Calculates residual for testing purposes....
! c   needs to run with double precision vectors to be useful.....
! c******************************************************************
!       subroutine testinv(Phi,resmax,itercg,am,imass,x,aden,ndiag)
!       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
!       parameter(kferm=4*kthird*kvol)
!       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
!       common/para/bbb,am3,ibound
!       complex Phi(kferm)
!       complex(dp) vtild(kferm)
!       complex(dp) x3(kferm)
!       complex(dp) x(kferm,ndiag),x1(kferm),x2(kferm)
!       complex(dp) u
! c     complex vtild(kferm)
! c     complex x3(kferm)
! c     complex x(kferm,ndiag),x1(kferm),x2(kferm)
! c      complex u
!       real(dp) residual
!       real(dp) aden(ndiag)
! c
!       write(6,111)
! 111   format(' Hi from testinv')
! c
!       resmax=0.0
! c
!       do idiag=1,ndiag
!       residual=0.0
!       do i=1,kferm
!       x3(i)=x(i,idiag)
!       enddo
!       call dslash(x2,x3,u,am,imass)
!       call dslashd(x1,x2,u,am,imass)
!       do i=1,kferm
!       vtild(i)=x1(i)+aden(idiag)*x3(i)-Phi(i)
!       residual=residual+conjg(vtild(i))*vtild(i)
!       enddo
! c     residual=sqrt(residual)
!       if(residual.gt.resmax) resmax=residual
! c
!       write(6,*) idiag, 'itercg = ',itercg, ' residual = ',residual
!       enddo
! c
!       resmax=sqrt(resmax)
! c
!       return
!       end

  subroutine sread
    use random
    use gauge
#ifdef MPI
    use comms
    integer :: mpi_fh
    integer :: status(mpi_status_size)
    integer :: ierr

    call MPI_File_Open(comm, 'con', MPI_Mode_Rdonly, &
         & MPI_Info_Null, mpi_fh, ierr)
! Get the configuration
    call MPI_File_Set_View(mpi_fh, 0_8, MPI_Real, mpiio_type, "native", &
         & MPI_Info_Null, ierr)
    call MPI_File_Read_All(mpi_fh, theta, 3*ksizex_l*ksizey_l*ksizet_l, &
         & MPI_Real, status, ierr)
    call MPI_File_Close(mpi_fh, ierr)
! Get the see,ierrd
    if (ip_global .eq. 0) then
      print *, "configuration file read."
      open (unit=10, file='con', status='old', form='unformatted', access='stream')
      !print*,"FSEEK CALL COMMENTED OUT, THIS WILL FAIL"
      call fseek(10, 3*ksize*ksize*ksizet*4 + 4, 0)
      read (10) seed
      close (10)
    end if
#else
    open (unit=10, file='con', status='old', form='unformatted')
    read (10) theta, seed
    close (10)
    print *, "configuration file read."
#endif
    return
  end subroutine sread
!
  subroutine swrite
    use random
    use gauge
#ifdef MPI
    use comms
    integer :: mpi_fh
    integer :: status(mpi_status_size)
    integer :: ierr

! Write theta
    call MPI_File_Open(comm, 'con', MPI_Mode_Wronly + MPI_Mode_Create, &
         & MPI_Info_Null, mpi_fh, ierr)
    call MPI_File_Set_View(mpi_fh, 0_8, MPI_Real, mpiio_type, "native", &
         & MPI_Info_Null, ierr)
    call MPI_File_Write_All(mpi_fh, theta, 3*ksizex_l*ksizey_l*ksizet_l, &
         & MPI_Real, status, ierr)
    call MPI_File_Close(mpi_fh, ierr)

! Write seed in serial
    if (ip_global .eq. 0) then
      open (unit=31, file='con', status='old', form='unformatted', access='stream')
      !print*,"FSEEK CALL COMMENTED OUT, THIS WILL FAIL"
      call fseek(31, 3*ksize*ksize*ksizet*4 + 4, 0)
! Manually compute the effective record length to be compatible with serial Fortran
      write (31) seed, 3*ksize*ksize*ksizet*4 + 8
      close (31)
      open (unit=40, file='random_seed', status='replace')
      write (40, *) seed
      close (40)

    end if
#else
    open (unit=31, file='con', status='unknown', form='unformatted')
    write (31) theta, seed
    close (31)
#endif
    return
  end subroutine swrite

  ! Tries to read seed from file
  subroutine readseed(globalseed)
#ifdef MPI
    use comms
#endif
    implicit none
    ! if the seed file is not correctly read, we do not want to change
    ! the vaule of globalseed passed as input
    real(dp), intent(inout) :: globalseed
    integer :: seedreadstatus
#ifdef MPI
    integer :: ierr
    if (ip_global .eq. 0) then
#endif
      open (unit=40, file='random_seed', iostat=seedreadstatus, status='old')
      ! if 'random_seed' can be opened, read from it.
      if (seedreadstatus .eq. 0) then
        read (40, *) globalseed
        close (40)
      endif
#ifdef MPI
    endif
    call MPI_Bcast(globalseed, 1, MPI_Double_Precision, 0, comm, ierr)
#endif
  end subroutine readseed

  ! Saves seed into file
  subroutine saveseed(seedtosave)
#ifdef MPI
    use comms
#endif
    implicit none
    real(dp), intent(in) :: seedtosave

#ifdef MPI
    if (ip_global .eq. 0) then
#endif
      open (unit=40, file='random_seed')
      write (40, *) seedtosave
      close (40)
#ifdef MPI
    endif
#endif
  end subroutine saveseed
!
  subroutine init(nc)
    use random
    use gauge
    use dirac
!*******************************************************************
!     sets initial values
!     nc=0 cold start
!     nc=1 hot start
!     nc<0 no initialization
!*******************************************************************
    integer, intent(in) :: nc
!     complex one,zi
    complex(dp), parameter :: one = (1.0, 0.0), zi = (0.0, 1.0)
    integer :: ix, iy, it, mu
    real :: g
!
!*******************************************************************
!  calculate constants
!*******************************************************************
!      call addrc
!*******************************************************************
!    setup Dirac algebra
!*******************************************************************
!
!     gamma_1
!
    gamval(1, 1) = -zi
    gamval(1, 2) = -zi
    gamval(1, 3) = zi
    gamval(1, 4) = zi
!
    gamin(1, 1) = 4
    gamin(1, 2) = 3
    gamin(1, 3) = 2
    gamin(1, 4) = 1
!
!     gamma_2
!
    gamval(2, 1) = -one
    gamval(2, 2) = one
    gamval(2, 3) = one
    gamval(2, 4) = -one
!
    gamin(2, 1) = 4
    gamin(2, 2) = 3
    gamin(2, 3) = 2
    gamin(2, 4) = 1
!
!     gamma_3
!
    gamval(3, 1) = -zi
    gamval(3, 2) = zi
    gamval(3, 3) = zi
    gamval(3, 4) = -zi
!
    gamin(3, 1) = 3
    gamin(3, 2) = 4
    gamin(3, 3) = 1
    gamin(3, 4) = 2
!
!     gamma_4
!
    gamval(4, 1) = one
    gamval(4, 2) = one
    gamval(4, 3) = -one
    gamval(4, 4) = -one
!
    gamin(4, 1) = 1
    gamin(4, 2) = 2
    gamin(4, 3) = 3
    gamin(4, 4) = 4
!
!     gamma_5 = gamma_1 * gamma_2 * gamma_3 * gamma_4
!
    gamval(5, 1) = -one
    gamval(5, 2) = -one
    gamval(5, 3) = -one
    gamval(5, 4) = -one
!
    gamin(5, 1) = 3
    gamin(5, 2) = 4
    gamin(5, 3) = 1
    gamin(5, 4) = 2
!
!     gamma_4 * gamma_5 (called gamma_3 gamma_5 in notes)
    gamval(6, 1) = -one
    gamval(6, 2) = -one
    gamval(6, 3) = one
    gamval(6, 4) = one
!
    gamin(6, 1) = 3
    gamin(6, 2) = 4
    gamin(6, 3) = 1
    gamin(6, 4) = 2
!
!
    gamval = gamval*akappa
!
    if (nc .lt. 0) return
!
!     initialize gauge fields
!
    if (nc .eq. 1) goto 40
!     (else cold start)
    theta = 0.0
    return
!
40  continue
    g = 0.05
    do mu = 1, 3
      do it = 1, ksizet_l
        do iy = 1, ksizey_l
          do ix = 1, ksizex_l
!               theta(ix, iy, it, mu) = 2.0 * g * rranf(ix, iy, it) - 1.0
            theta(ix, iy, it, mu) = 2.0*g*rano(yran, idum, ix, iy, it) - 1.0
          enddo
        enddo
      enddo
    enddo
    return
  end subroutine init
!******************************************************************
!   calculate compact links from non-compact links
!******************************************************************
  subroutine coef(u, theta)
    use comms4, only: start_halo_update_4
    use comms_common, only: ip_t
    use comms, only: complete_halo_update
    implicit none
!
    complex(dp), intent(out) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: theta(ksizex_l, ksizey_l, ksizet_l, 3)
#ifdef MPI
    integer, dimension(12) :: reqs_u
#endif

!     u(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = exp(cmplx(0.0, theta))
    u(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = (1.0 + cmplx(0.0, theta))
!
!  anti-p.b.c. in timelike direction
    if (ibound .eq. -1 .and. ip_t .eq. (np_t - 1)) then
      u(:, :, ksizet_l, 3) = -u(:, :, ksizet_l, 3)
    end if
!
!!!!    call complete_halo_update_4(3, u)
#ifdef MPI
    call start_halo_update_4(3, u, 3, reqs_u)
    call complete_halo_update(reqs_u)
#else
    call update_halo_4(3, u)
#endif
    return
  end subroutine coef
end module dwf3d_lib
