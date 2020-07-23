program full_md
  use gdbhook
  use random
  use gaussian
  use remez
  use remezg
  use trial, ut => u, thetat => theta
  use gauge
  use vector, X1 => X
  use gforce
  use counters
  use dum1
  use comms
  use measure_module
  use qmrherm_module, only: qmrherm, qmrhprint => printall
  use timer, only: timeinit => initialise, get_time_from_start
  use evolution, only: evolve_theta_pp, initialise_phi_1flavour, initialise_pp
  use dwf3d_lib
  implicit none
  !real, parameter :: respbp=1.0e-6, rescgg=1.0e-6
  !real, parameter :: rescga=1e-9
  !real, parameter :: rescgm=1e-9
  double precision :: t1i, t2i
  !integer, parameter :: itermax=1000
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)!
  complex(dp) :: Xresult(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real(dp) :: H0, H1, S0, S1, dH, dS, hg, hp
  real :: action, paction, gaction
  real :: vel2, x, atraj
  real :: dt, am, y, traj, proby
  integer :: imass, iter, iterl, iter2, i, ia, idirac, ithird
  integer :: isweep, itercg, mu
  !*******************************************************************
  !     variables to keep track of MPI requests
  !*******************************************************************
#ifdef MPI
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
  call gdb_wait()
  call timeinit
#endif
  !*******************************************************************
  !     end of input
  !*******************************************************************
  !*******************************************************************
  !     check qmrherm is going to be OK
  !*******************************************************************
  if (ip_global .eq. 0) then
    open (unit=7, file='output', status='unknown')
    open (unit=98, file='control', status='unknown')
  end if

  dt = 0.12 ! lenght of md step
  beta = 0.4
  am3 = 1.0
  am = 0.05
  imass = 3
  iterl = 4 ! nsteps of md trajectory
  iter2 = 1

  ! set a new seed by hand...

  seed = 4139764973254.0 ! this should be a double precision number
  ! ending with d0, but we do not care
  call init_random(seed)

  call init(1) ! hot random start
  !  mock Remez coefficients
  anum2(0) = 71
  anum4(0) = 71
  anum2g(0) = 71
  anum4g(0) = 71
  do i = 1, ndiag
    anum2(i) = -2.0879703342431707e-08
    aden2(i) = 1.0032174381543456e-05
    anum4(i) = -2.0879703342431707e-08
    aden4(i) = 1.0032174381543456e-05
  enddo
  do i = 1, ndiagg
    anum2g(i) = -2.0879703342431707e-08
    aden2g(i) = 1.0032174381543456e-05
    anum4g(i) = -2.0879703342431707e-08
    aden4g(i) = 1.0032174381543456e-05
  enddo
  bnum2 = anum2
  bnum4 = anum4
  bnum2g = anum2g
  bnum4g = anum4g

  !*******************************************************************
  !     print heading
  !*******************************************************************
  traj = iterl*dt
  proby = 1.0/float(iterl)
  if (ip_global .eq. 0) then
    write (*, *) "Benchmark - fake rational approximation data"
    write (7, 9001) ksize, ksizet, kthird, Nf, dt, traj, ndiag, ndiagg, &
      & iter2, beta, am3, am, imass
9001 format(' ksize=', i3, ' ksizet=', i3, / &
           , ' kthird=', i3, / &
           , ' Nf =', i3, / &
           , ' time step: dt=', f6.4, ' trajectory length=', f9.6, / &
           , ' Remez ndiag: action =', i3, ' guidance=', i3, / &
           , ' # trajectories=', i6, ' beta=', f9.6, / &
           , ' am3=', f6.4, ' am=', f6.4/ &
           , ' imass=', i2)
#ifdef MPI
    write (7, 9002) NP_X, NP_Y, NP_T, ksizex_l, ksizey_l, ksizet_l
9002 format(" NP_X=", i3, " NP_Y=", i3, " NP_T=", i3, / &
& " ksizex_l=", i3, " ksizey_l=", i3, " ksizet_l=", i3)
#endif
    !     write(6,9004) rescgg,rescga,respbp
    write (7, 9004) rescgg, rescga, respbp
9004 format(' Stopping residuals: guidance: ', e11.4, ' acceptance: ', &
&     e11.4, ' estimator: ', e11.4)
    !     write(6,9044) rescgm
    write (7, 9044) rescgm
9044 format(' Stopping residuals: meson: ', e11.4)
    call rranget(seed, 1, 1, 1)
    ! c     write(6,*) 'seed: ', seed
    write (7, *) 'seed: ', seed
  end if
  !*******************************************************************
  !       initialize for averages
  !*******************************************************************

  call init_counters()

  !*******************************************************************
  !     start of classical evolution
  !*******************************************************************
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

    !! >> initialise_phi
    call initialise_phi_1flavour(phi, am, imass)

    !! >> initialise_pp
    call initialise_pp(pp)
    ! Start computing time (including hamiltonian calls)

    call hamilton(Phi, H0, hg, hp, S0, rescga, am, imass)
    if (isweep .eq. 1) then
      action = real(S0)/kvol
      gaction = real(hg)/kvol
      paction = real(hp)/kvol
    endif

#ifdef MPI
    t1i = MPI_Wtime()
#endif
    call evolve_theta_pp(iterl, dt, pp, phi, rescgg, am, imass, itot)

    !**********************************************************************
    !  Monte Carlo step: accept new fields with probability=
    !              min(1,exp(H0-H1))
    !**********************************************************************
501 continue
#ifdef MPI
    t2i = MPI_Wtime()
    if (ip_global .eq. 0) then
#endif
      write (6, *) "End of MD."
      time_measurement: block
        double precision :: dt
        dt = (t2i - t1i)/iterl
        print *, "Time per MD iteration:", dt
      end block time_measurement
      print *, dt
#ifdef MPI
    endif
#endif
    call hamilton(Phi, H1, hg, hp, S1, rescga, am, imass)
    dH = H0 - H1
    dS = S0 - S1
    if (ip_global .eq. 0) then
      write (98, *) dH, dS
      write (6, *) "dH,dS ", dH, dS
    end if
    y = exp(real(dH))
    y_average = y_average + y
    ysq_average = ysq_average + y*y

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
600 continue
    if (ip_global .eq. 0) then
      write (11, *) isweep, gaction, paction
    end if
    action_average = action_average + action
    vel2 = sum(pp*pp)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, vel2, 1, MPI_Real, MPI_Sum, comm, ierr)
#endif
    vel2 = vel2/(3*kvol)
    vel2a = vel2a + vel2
    !
    !     uncomment to disable measurements
    !     goto 601
    !666    continue
    if ((isweep/iprint)*iprint .eq. isweep) then
      thetat = theta
      call coef(ut, thetat)
      call measure(pbp, respbp, ancgm, am, imass)
      !         call meson(rescgm,itercg,ancgm,am,imass)
      pbp_average = pbp_average + pbp
      ancgm_average = ancgm_average + ancgm
      ipbp = ipbp + 1
      !        write(11,*) pbp
      if (ip_global .eq. 0) then
        write (6, *) isweep, 'pbp:', pbp, ancgm
      endif
    endif

  end do
  !*******************************************************************
  !     end of main loop
  !*******************************************************************
  call final_averages(Nf, iter2)
  atraj = dt*itot/iter2
  !*******************************************************************
  !     print global averages
  !*******************************************************************
  if (ip_global .eq. 0) then
    write (7, 9022) iter2, naccp, atraj, y_average, ysq_average, &
      & ancg, ancgpv, ancgh, ancghpv, ancgf, ancgfpv, ancgpf, ancgpfpv, &
      & pbp_average, ancgm_average, vel2a, action_average
9022 format(' averages for last ', i6, ' trajectories', /  &
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
9024 format(1x)
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
end program full_md

