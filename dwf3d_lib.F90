module dwf3d_lib
  use params
  use mpi
  implicit none
  save

  ! Random numbers
  real(dp) :: seed

! Useful constants
  real, parameter :: One = 1.0
contains

  pure subroutine verify_kernel_choice()
    ! Measurement flag
#if defined(MEASURE_SHAMIR) && defined(MEASURE_WILSON)
    Error: Must specify only one of MEASURE_SHAMIR or MEASURE_WILSON
#endif
#if !defined(MEASURE_SHAMIR) && !defined(MEASURE_WILSON)
    Error: Must specify one of MEASURE_SHAMIR or MEASURE_WILSON
#endif
    
    ! Production flag
#if defined(GENERATE_WITH_SHAMIR) && defined(GENERATE_WITH_WILSON)
    Error: Must specify only one of GENERATE_WITH_SHAMIR or GENERATE_WITH_WILSON
#endif
#if !defined(GENERATE_WITH_SHAMIR) && !defined(GENERATE_WITH_WILSON)
    Error: Must specify one of GENERATE_WITH_SHAMIR or GENERATE_WITH_WILSON
#endif
  end subroutine

  subroutine dwf3d_main
    use gdbhook
    use random
    use gaussian
    use remez
    use remezg
    use remez_common_subroutines
    use trial, ut => u, thetat => theta
    use gauge, only: theta, coef
    use vector, X1 => X
    use counters ! ALL
    use dum1
    use comms
    use measure_module
    use qmrherm_module, only: qmrherm, qmrhprint => printall
    use timer, only: timeinit => initialise, get_time_from_start
    use evolution, only: evolve_theta_pp, initialise_phi_1flavour, initialise_pp
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
!    action_average: running average of total action
!
!                                               SJH February 2017

!    The code now runs up to a fixed numner iter2_read of trajectories, but run time is now the limiting factor (see
!    walltime_seconds variable.
!*******************************************************************
    implicit none
    complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)!
    real(dp) :: H0, H1, S0, S1, dH, dS, hg, hp
    real :: action, paction, gaction
    real :: vel2, x, atraj
    real :: dt, am, y, traj
    integer :: imass, iterl, iter2, iter2_read
    integer :: walltimesec
    logical :: program_status_file_exists
    
#ifdef MPI
!     variables to keep track of MPI requests
    integer :: ierr
    integer :: itercg
    ! real :: sumvalue, maxvalue

#endif
    call verify_kernel_choice()
    ibound = -1
    qmrhprint = .true.
#ifdef MPI
    call init_MPI
    call gdb_wait()
    call timeinit
#endif

    if (ip_global .eq. 0) then
      open (unit=7, file='output', status='unknown', action='write', &
            position='append')
      open (unit=98, file='control', status='unknown', action='write', &
            position='append')
    end if

    open (unit=25, file='midout', status='old', action='read')

    read (25, *) dt, beta, am3, am, imass, iterl, iter2_read, walltimesec
    close (25)

    ! verify imass value
    if ((imass .ne. 1) .and. (imass .ne. 3) .and. (imass .ne. 5)) then
      print *, 'ERROR: imass must be one of 1, 3 or 5'
#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#endif
      call exit(1)
    else if (imass .eq. 5)
      print *, 'WARNING: imass of 5 may not be supported.'
    end if

! set a new seed by hand...

    block
      logical :: success
      call readseed(seed, success)
      if (.not. success) then
#ifdef MPI
        if (ip_global .eq. 0) then
#endif
          print *, 'Reading seed file failed, using default value'
#ifdef MPI
        endif
#endif
        seed = 4139764973254.0 ! this should be a double precision number
        ! ending with d0, but we do not care
      endif
    end block
    action = -1

    if (ip_global .eq. 0) then
      write (7, *) 'seed: ', seed
    end if

#ifdef MPI
    call seed_rank(seed)
#endif

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

    call read_remez_file('remez2', ndiag, anum2, bnum2, aden2, bden2)
    call read_remez_file('remez2g', ndiagg, anum2g, bnum2g, aden2g, bden2g)

    call read_remez_file('remez4', ndiag, anum4, bnum4, aden4, bden4)
    call read_remez_file('remez4g', ndiagg, anum4g, bnum4g, aden4g, bden4g)

!*******************************************************************
!     print heading
!*******************************************************************
    traj = iterl*dt
!*******************************************************************
!       initialize for averages
!*******************************************************************

    call init_counters()

!*******************************************************************
!     start of classical evolution
!*******************************************************************
    classical_evolution: block

      real :: run_time, time_per_md_step ! conservative estimates
      real :: measurement_time, total_md_time
      integer :: isweep, isweep_total_start

      measurement_time = 100 ! an arbitrary value (that should be conservative)
      total_md_time = 0
      inquire (file='program_status', exist=program_status_file_exists)
      if (program_status_file_exists) then
        open (unit=53, file='program_status', status='old', action='read')
        read (53, *) isweep_total_start, measurement_time
        close (53)
      else
        isweep_total_start = 0
      endif

      do isweep = 1, iter2_read

#ifdef MPI
        if (ip_global .eq. 0) then
#endif
          write (6, "(A7,I6,A4,I6,A12,I6,A1)") 'Isweep ', isweep, ' of ', iter2_read, ' (start from', isweep_total_start, ')'
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
        run_time = get_time_from_start()

        total_md_time = total_md_time - run_time
        call hamilton(Phi, H0, hg, hp, S0, rescga, am, imass)

        if (isweep .eq. 1) then
          action = real(S0)/kvol
          gaction = real(hg)/kvol
          paction = real(hp)/kvol
        endif

        call evolve_theta_pp(iterl, dt, pp, phi, rescgg, am, imass, itot)

!**********************************************************************
!  Monte Carlo step: accept new fields with probability=
!              min(1,exp(H0-H1))
!**********************************************************************

        call hamilton(Phi, H1, hg, hp, S1, rescga, am, imass)
        dH = H0 - H1
        dS = S0 - S1
        if (ip_global .eq. 0) then
          write (98, *) isweep_total_start, isweep, dH, dS
          write (6, *) "dH,dS ", dH, dS
        end if
        y = exp(real(dH))
        y_average = y_average + y
        ysq_average = ysq_average + y*y

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
          open (unit=11, file='fort.11', action='write', position='append')
          write (11, *) isweep_total_start, isweep, gaction, paction
          close (11)
        end if
        if (action .eq. -1) then
          print *, "action was not correctly initialised."
          stop
        end if
        action_average = action_average + action
        vel2 = sum(pp*pp)
#ifdef MPI
        call MPI_AllReduce(MPI_In_Place, vel2, 1, MPI_Real, MPI_Sum, comm, ierr)
        vel2 = vel2/np_third
#endif
        vel2 = vel2/(3*kvol)
        vel2a = vel2a + vel2

        ! Including also hamiltonian call time
        run_time = get_time_from_start()

        total_md_time = total_md_time + run_time
        time_per_md_step = total_md_time/itot

!     uncomment to disable measurements
!     goto 601
!666    continue

        if (mod((isweep + isweep_total_start), iprint) .eq. 0) then
          thetat = theta
          call coef(ut, thetat)
!          call measure(pbp, respbp, ancgm, am, imass, isweep + isweep_total_start)
          call meson(rescgm,itercg,ancgm,am,imass, isweep + isweep_total_start)
          call fermion(rescgm,itercg,ancgm,am,imass, isweep + isweep_total_start)
          pbp_average = pbp_average + pbp
          ancgm_average = ancgm_average + ancgm
          ipbp = ipbp + 1
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            write (6, *) isweep, 'pbp:', pbp, ancgm
#ifdef MPI
          endif
#endif
          measurement_time = get_time_from_start()
          measurement_time = measurement_time - run_time
        endif
!
        if (mod((isweep + isweep_total_start), icheckpoint) .eq. 0) then
          call rranget(seed, 1, 1, 1)
          if (iwrite .eq. 1) then
            call swrite(isweep + isweep_total_start)
          endif
          !  now dealt with properly in measure_module
          !flush (100)
          !flush (200)
        endif

        keep_running_check: block
          real :: time_for_next_iteration
          real :: time_left
          logical :: stop_file_exists

          time_for_next_iteration = time_per_md_step*4*iterl*2

          if (mod((isweep + isweep_total_start + 1), iprint) .eq. 0) then
            time_for_next_iteration = time_for_next_iteration + measurement_time
          endif

          run_time = get_time_from_start()

#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            print *, 'Expected next run time with safety margin:', &
              run_time + time_for_next_iteration, ' of ', walltimesec
            print *, 'Run so far for: ', run_time
#ifdef MPI
          endif
#endif
          time_left = walltimesec - (run_time + time_for_next_iteration)
          call MPI_Bcast(time_left,1,MPI_Real,0,comm, ierr)

          if (time_left .lt. 0) then
#ifdef MPI
            if (ip_global .eq. 0) then
#endif
              print *, 'Time left insufficient, quitting.'
              print *, 'Run:', isweep, ' started:', isweep_total_start
#ifdef MPI
            endif
#endif
            exit
          endif

          inquire (file='stop', exist=stop_file_exists)
          if (stop_file_exists) then
#ifdef MPI
            if (ip_global .eq. 0) then
#endif
              print *, 'Found "stop" file: stopping now.'
              print *, 'Run:', isweep, ' started:', isweep_total_start
#ifdef MPI
            endif
#endif
            exit
          endif

        end block keep_running_check
      end do
      iter2 = min(isweep, iter2_read)
      if (ip_global .eq. 0) then
        open (unit=53, file='program_status', status='unknown', action='write')
        write (53, *) iter2 + isweep_total_start, measurement_time
        close (53)
      endif
    end block classical_evolution
    !*******************************************************************
!     end of main loop
!*******************************************************************
    call final_averages(Nf, iter2)

    atraj = dt*itot/iter2
!*******************************************************************
!     print global averages
!*******************************************************************

    if (ip_global .eq. 0) then
      write (7, 9001) ksize, ksizet, kthird, Nf, dt, traj, ndiag, ndiagg, &
        iter2, iter2_read, beta, am3, am, imass
9001  format(' ksize=', i3, ' ksizet=', i3, / &
             , ' kthird=', i3, / &
             , ' Nf =', i3, / &
             , ' time step: dt=', f6.4, ' expected avg trajectory length=', f9.6, / &
             , ' Remez ndiag: action =', i3, ' guidance=', i3, / &
             , ' # trajectories=', i6, ' out of max', i6, ' beta=', f9.6, / &
             , ' am3=', f6.4, ' am=', f6.4/ &
             , ' imass=', i2)
#ifdef MPI
      write (7, 9002) NP_X, NP_Y, NP_T, NP_THIRD, ksizex_l, ksizey_l, ksizet_l, kthird_l
9002  format(" NP_X=", i3, " NP_Y=", i3, " NP_T=", i3, " NP_THIRD=", i3, / &
             " ksizex_l=", i3, " ksizey_l=", i3, " ksizet_l=", i3, " kthird_l=", i3)
#endif
!     write(6,9004) rescgg,rescga,respbp
      write (7, 9004) rescgg, rescga, respbp
9004  format(' Stopping residuals: guidance: ', e11.4, ' acceptance: ', &
           &     e11.4, ' estimator: ', e11.4)
!     write(6,9044) rescgm
      write (7, 9044) rescgm
9044  format(' Stopping residuals: meson: ', e11.4)
      call rranget(seed, 1, 1, 1)
      write (7, *) 'seed: ', seed
      write (7, 9022) iter2, naccp, atraj, y_average, ysq_average, &
           & ancg, ancgpv, ancgh, ancghpv, ancgf, ancgfpv, ancgpf, ancgpfpv, &
           & pbp_average, ancgm_average, vel2a, action_average
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
    endif

#ifdef MPI
    call MPI_Finalize(ierr)
#endif
!
  end subroutine dwf3d_main
!******************************************************************
!   Evaluation of Hamiltonian function
!******************************************************************
  subroutine hamilton(Phi, h, hg, hp, s, res2, am, imass, abort_on_max_reached_in)
    use remez
    use trial, only: theta, pp
    use dum1
    use counters, only: ancghpv, ancgh
    use comms
    use qmrherm_module, only: qmrherm
    implicit none
    complex(dp), intent(in) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, Nf)
    real(dp), intent(out) :: h, hg, hp, s
    real, intent(in) :: res2, am
    integer, intent(in) :: imass
    logical, intent(in), optional :: abort_on_max_reached_in
    complex(dp) :: Xresult(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real(dp) :: hf
    integer :: itercg, ia
    logical :: abort_on_max_reached
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
    hp = hp/np_third
#endif

    hg = 0.5*Nf*beta*sum(theta**2)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, hg, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
    hg = hg/np_third
#endif
    h = hg + hp

    if (present(abort_on_max_reached_in)) then
      abort_on_max_reached = abort_on_max_reached_in
    else
      abort_on_max_reached = .true.
    endif

! uncomment these lines to quench the fermions!
!     return
!
!  pseudofermion action is
!   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
!
    do ia = 1, Nf
      R = Phi(:, :, :, :, :, ia)

      call qmrherm(R, Xresult, res2, itercg, One, 1, anum4, aden4, ndiag, 0)
      ancghpv = ancghpv + float(itercg)
      call check_qmr_iterations(niterations=itercg, abort_on_max_reached=abort_on_max_reached)

      R = Xresult

      call qmrherm(R, Xresult, res2, itercg, am, imass, bnum2, bden2, ndiag, 0)
      ancgh = ancgh + float(itercg)
      call check_qmr_iterations(niterations=itercg, abort_on_max_reached=abort_on_max_reached)

      hf = hf + sum(real(conjg(R(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) &
                         *Xresult(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))
    enddo

#ifdef MPI
! hf is built up from zero during the loop so only needs to be summed across
! all partitions at this point
    call MPI_AllReduce(MPI_In_Place, hf, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
#endif

    h = hg + hp + hf
    s = hg + hf
!
    return
  end subroutine hamilton

  subroutine sread(success_out)
    use random
    use gauge
#ifdef MPI
    use comms
    implicit none
    logical, optional, intent(out) :: success_out
    logical :: success
    integer :: mpi_fh
    integer :: status(mpi_status_size)
    integer :: ierr

    inquire (file='con', exist=success)

    if (success) then

      call MPI_File_Open(comm, 'con', MPI_Mode_Rdonly, &
                         MPI_Info_Null, mpi_fh, ierr)
! Get the configuration
      call MPI_File_Set_View(mpi_fh, 0_8, MPI_Real, mpiio_type, "native", &
                             MPI_Info_Null, ierr)
      call MPI_File_Read_All(mpi_fh, theta, 3*ksizex_l*ksizey_l*ksizet_l, &
                             MPI_Real, status, ierr)
      call MPI_File_Close(mpi_fh, ierr)
! Get the see,ierrd
      if (ip_global .eq. 0) then
        open (unit=10, file='con', status='old', form='unformatted', access='stream')
        !print*,"FSEEK CALL COMMENTED OUT, THIS WILL FAIL"
        call fseek(10, 3*ksize*ksize*ksizet*4 + 4, 0)
        read (10) seed
        close (10)
        print *, "configuration file read."
      end if
#else
      open (unit=10, file='con', status='old', form='unformatted')
      read (10) theta, seed
      close (10)
      print *, "configuration file read."
#endif
    else
      if ((.not. success) .and. (.not. present(success_out))) then
        print *, "Configuration file 'con' not found!"
        stop
      endif
    endif
    if (present(success_out)) then
      success_out = success
    endif

  end subroutine sread

  subroutine swrite(traj_id)
    use random
    use gauge
#ifdef MPI
    use comms
    implicit none
    integer, intent(in), optional :: traj_id

    integer :: mpi_fh
    integer :: status(mpi_status_size)
    integer :: ierr
    character(len=20) :: configuration_filename

    if (present(traj_id)) then
      write (configuration_filename, '(A4,I5.5)') 'con.', traj_id
    else
      write (configuration_filename, '(A3)') 'con'
    endif
! Write theta
    call MPI_File_Open(comm, trim(configuration_filename), &
                       MPI_Mode_Wronly + MPI_Mode_Create, MPI_Info_Null, mpi_fh, ierr)
    call MPI_File_Set_View(mpi_fh, 0_8, MPI_Real, mpiio_type, "native", &
                           MPI_Info_Null, ierr)
    call MPI_File_Write_All(mpi_fh, theta, 3*ksizex_l*ksizey_l*ksizet_l, &
                            MPI_Real, status, ierr)
    call MPI_File_Close(mpi_fh, ierr)

! Write seed in serial
    if (ip_global .eq. 0) then
      open (unit=31, file=trim(configuration_filename), status='old', &
            form='unformatted', access='stream')
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
    open (unit=31, file=trim(configuration_filename), status='unknown', &
          form='unformatted')
    write (31) theta, seed
    close (31)
#endif
    return
  end subroutine swrite

  ! Tries to read seed from file
  subroutine readseed(globalseed, success)
#ifdef MPI
    use comms
#endif
    implicit none
    ! if the seed file is not correctly read, we do not want to change
    ! the vaule of globalseed passed as input
    real(dp), intent(inout) :: globalseed
    logical, intent(out) :: success
    integer :: seedreadstatus
#ifdef MPI
    integer :: ierr
    if (ip_global .eq. 0) then
#endif
      open (unit=40, file='random_seed', iostat=seedreadstatus, status='old')
      ! if 'random_seed' can be opened, read from it.
      if (seedreadstatus .eq. 0) then
        read (40, *) globalseed
        success = .true.
        close (40)
      else
        success = .false.
      endif
#ifdef MPI
    endif
    call MPI_Bcast(globalseed, 1, MPI_Double_Precision, 0, comm, ierr)
    call MPI_Bcast(success, 1, MPI_Logical, 0, comm, ierr)
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

#ifdef MPI
! This subroutine takes the global seed and changes it so that it is
! different for each rank.
! Makes also sure that is positive and greater than 1.

  subroutine seed_rank(seed)
    use comms_common, only: ip_global
    real(dp), intent(inout) :: seed

#ifdef SITERANDOM
    if (NP_THIRD .ne. 1) then
      print *, "SITERANDOM is not fully compatible with NP_THIRD > 1"
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      stop
    endif
#endif
    seed = seed + ip_global

    if (seed .lt. 1) then
      seed = abs(seed) + 1
    endif
  end subroutine seed_rank
#endif

  subroutine init(nc)
    use gammamatrices, only: init_gammas
    use gauge
    use random
!*******************************************************************
!     sets initial values
!     nc=0 cold start
!     nc=1 hot start
!     nc<0 no initialization
!*******************************************************************
    implicit none
    integer, intent(in) :: nc
    integer :: nc_temp
    integer :: ix, iy, it, mu
    real :: g

    call init_gammas()
    block
      logical :: success
      if ((istart .lt. 0) .or. (iread .eq. 1)) then
        call sread(success)
        if (.not. success) then
#ifdef MPI
          if (ip_global .eq. 0) then
#endif
            print *, 'Reading con file failed, starting from random conf'
#ifdef MPI
          endif
#endif
          nc_temp = 1
        else
          return
        endif
      else
        nc_temp = nc
      endif
    end block

#ifdef MPI
    if (ip_global .eq. 0) then
#endif
      print *, 'Initialising gauge conf, nc: ', nc_temp
#ifdef MPI
    endif
#endif

    select case (nc_temp)
    case (-1)
      return
    case (0)
      theta = 0.0
      return
    case (1)
      g = 0.05
      do mu = 1, 3
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              theta(ix, iy, it, mu) = 2.0*g*rano(yran, idum, ix, iy, it) - 1.0
            enddo
          enddo
        enddo
      enddo
      return
    end select
  end subroutine init

  subroutine check_qmr_iterations(niterations, abort_on_max_reached)
    use params, only: max_qmr_iters
#ifdef MPI
    use mpi
    integer :: ierr
#endif
    integer, intent(in) :: niterations
    logical, intent(in) :: abort_on_max_reached

    if (niterations .eq. max_qmr_iters) then
      if (abort_on_max_reached) then
        print *, "ERROR: Max QMR iterations reached."
#ifdef MPI
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#else
        stop
#endif
      else
        print *, "WARNING: Max QMR iterations reached - NOT Aborting, as requested."
      endif
    endif
  end subroutine check_qmr_iterations

end module dwf3d_lib

