program dry_convection
    ! Import modules
    ! use constants_module
    use finitedifference
    use namelist_utilities
    use csv_writer
    use matrix_utilities
    use boundaries_ex5
    use T_inits
    use poisson_solver
    use poisson_duetsch
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx=4, ny=4                    ! number of grid points in x and y
    integer :: io_error                      ! error code for file IO
    integer :: k                            ! iteration counter
    integer :: nsteps                        ! number of time steps
    real :: a_adv, a_diff                    ! advection and diffusion coefficients
    real :: total_time                       ! total simulation time
    real :: max_err                          ! maximum error
    real :: Ra                               ! Rayleigh number
    real(8) :: h                                ! grid spacing
    real(8) :: alpha = 1.0                           ! relaxation parameter
    REAL(8) :: res_rms                          ! root mean square residue 
    REAL(8) :: time = 0.0                     ! time. since h is change wee neet to track tiem
    character(len=100) :: T_ini_type        ! initial temperature profile type
    character(len=100) :: outputfilename    ! output file name
    real(8), allocatable:: T(:,:)               ! temperature profile
    real(8), allocatable:: T_new(:,:)           ! new temperature profile
    real(8), allocatable:: w(:,:)               ! vorticity
    real(8), allocatable:: psi(:,:)             ! stream function
    real(8), allocatable:: Tdx(:,:)             ! temperature derivative in x
    real(8), allocatable:: u(:,:)               ! x velocity
    real(8), allocatable:: v(:,:)               ! y velocity
    real(8), allocatable:: d2T(:,:)             ! second derivative(x,y) of temperature
    real(8), allocatable:: advection(:,:)       ! advection

    ! test derivetives:
    real(8), allocatable:: d2psi(:,:)           
    real(8), allocatable:: d2w(:,:)
    REAL(8) :: res_rms_duetsch                          ! root mean square residue 
    real(8), allocatable:: w_duetsch(:,:)               ! vorticity


    

    !read input parameters
    call read_namelist_ex7('data/namelist/ex_7.nml', nx, ny, a_adv, a_diff, total_time, max_err, Ra, T_ini_type)
    
    ! print namelist values
    print *, 'nx = ', nx
    print *, 'ny = ', ny
    print *, 'a_adv = ', a_adv
    print *, 'a_diff = ', a_diff
    print *, 'total_time = ', total_time
    print *, 'max_err = ', max_err
    print *, 'Ra = ', Ra
    print *, 'T_ini_type = ', T_ini_type


    ! Init Temperature arrays
    allocate(T(nx, ny), T_new(nx, ny), w(nx, ny), psi(nx, ny), Tdx(nx, ny), u(nx, ny), v(nx, ny))
    allocate(advection(nx, ny), d2T(nx, ny))

    ! Initialize temperature array
    T = 0.0
    T_new = 0.0
    w = 0.0
    psi = 0.0
    Tdx = 0.0
    u = 0.0
    v = 0.0
    d2T = 0.0
    advection = 0.0

    ! test derivetives:
    allocate(d2psi(nx, ny), d2w(nx, ny), w_duetsch(nx, ny))
    d2psi = 0.0
    d2w = 0.0
    w_duetsch = 0.0



    h = 1.0 / (ny - 1)

    nsteps = int(total_time / h)

    !call print_matrix(T)
    
    if (trim(T_ini_type) == 'rand') then
        ! Fill T with random numbers
        call RANDOM_NUMBER(T)
    else if (trim(T_ini_type) == 'cosine') then   !TODO:: change to cosine
        call initialize_temperature_cosine(T, nx, ny)
    else
        ! Handle unknown initialization type
        print *, 'Invalid initialization state:', trim(T_ini_type)
        stop
    end if


    ! save T for plotting
    open(unit=10, file=trim('data/ex_7/T.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    call write_to_csv_real8('data/ex_7/T.csv', T)
    close(10)


    ! save Tdx for plotting
    open(unit=10, file=trim('data/ex_7/Tdx.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    ! save w for plotting
    open(unit=10, file=trim('data/ex_7/w.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    ! save w for plotting
    open(unit=10, file=trim('data/ex_7/w.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    ! save psi for plotting
    open(unit=10, file=trim('data/ex_7/psi.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    ! save u for plotting
    open(unit=10, file=trim('data/ex_7/u.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    ! save psvi for plotting
    open(unit=10, file=trim('data/ex_7/v.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    ! save T_new for plotting
    open(unit=10, file=trim('data/ex_7/T_new.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)


    ! Save test derivertivs:
    
    ! save T_new for plotting
    open(unit=10, file=trim('data/ex_7/d2psi.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    
    ! save T_new for plotting
    open(unit=10, file=trim('data/ex_7/d2w.csv'), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    


    ! Integration loop
    do k = 1, int(nsteps)

        ! calc first derivertiv of temperature in x direction:
        !Compute dT/dx
        call dx(T, h, Tdx)

        call write_to_csv_real8('data/ex_7/Tdx.csv', Tdx)
        ! print max element
        print *, 'max Tdx = ', MAXVAL(ABS(Tdx))


        ! Determine 𝜔 from Ra dT/dx using the Poisson solver (Equation 2)
        res_rms = Vcycle_2DPoisson(w, -Ra*Tdx ,h, alpha) ! TODO:: check if - is correct
        print *, 'res_rms verena of w = ', res_rms
        call write_to_csv_real8('data/ex_7/w.csv', w)

        !Duetsch Sol:
        res_rms_duetsch = Vcycle_2DPoisson_duetsch(w_duetsch, -Ra*Tdx ,h, alpha) ! TODO:: check if - is correct
        print *, 'res_rms duetsch of w = ', res_rms_duetsch


        call FiniteDifference2D_real8(w/(-Ra), h, d2w)
        call write_to_csv_real8('data/ex_7/d2w.csv', d2w)


        ! Compute stream function 𝜓 from 𝜔 using the Poisson solver (Equation 3)
        res_rms = Vcycle_2DPoisson(psi, -w, h, alpha)  ! TODO:: check if - is correct
        call write_to_csv_real8('data/ex_7/psi.csv', psi)

        call FiniteDifference2D_real8(-psi, h, d2psi)
        call write_to_csv_real8('data/ex_7/d2psi.csv', d2psi)



        ! Compute the wind speeds 𝑢 and 𝑣 from 𝜓
        call dy(psi, h, u)
        call write_to_csv_real8('data/ex_7/u.csv', u)
        call dx(psi, h, v)
        call write_to_csv_real8('data/ex_7/v.csv', v)

        ! Compute the time step from a_adv, a_diff and the maximum wind speed in the
        ! we have no kappa this time like in ex5
        h = MIN(a_diff*h**2, a_adv*h/MAX(MAXVAL(ABS(u)), MAXVAL(ABS(v))))
        ! TODO: Müsste sich hier nciht auch nsteps änder?


        ! Compute ∇2𝑇 and 𝑣⃑ * ∇𝑇 using the subroutines from Exercise 5
        ! todo WHAT ARE DIM OF  𝑣⃑ * ∇𝑇 warum zweite deriv ein dim?

        ! ∇2𝑇:
        call FiniteDifference2D_real8(T, h, d2T)
        
        ! 𝑣⃑ * ∇𝑇 = advection
        call advection2D_real8(T, h, u, v, advection)



        ! calc next T
        !T_new = T + h* (d2T + advection)

        T = T_new
        call write_to_csv_real8('data/ex_7/T.csv', T)



        ! track time. not relevant for calc but since h is changing we need to track time
        time = time + h
        













        ! Maxiter to protect from to long runtimes and mage num instabel behavior debuggin easier
        if (k > 1) then
            exit
        end if
    end do



    







    print *, 'I am a zebra'





    
end program dry_convection











