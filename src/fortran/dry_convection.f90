program dry_convection
    ! Import modules
    ! use constants_module
    use finitedifference
    use namelist_utilities
    use csv_writer
    use matrix_utilities
    use boundaries_ex7
    use T_inits
    !use poisson_solver
    use poisson_solver_utilities
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx=4, ny=4                    ! number of grid points in x and y
    integer :: io_error                      ! error code for file IO
    integer :: k                            ! iteration counter
    integer :: nsteps                        ! number of time steps
    integer :: i, j                          ! loop indices
    real :: a_adv, a_diff                    ! advection and diffusion coefficients
    real :: total_time                       ! total simulation time
    real :: max_err                          ! maximum error
    real :: Ra                               ! Rayleigh number
    real(8) :: h                             ! grid spacing
    real(8) :: dt                            ! time step
    real(8) :: alpha = 1.                         ! relaxation parameter
    REAL(8) :: res_rms                          ! root mean square residue 
    REAL(8) :: time = 0.0                     ! time. since h is change wee neet to track time
    REAL(8) :: f_norm                          ! norm used for convergence check of poisson solver
    character(len=100) :: T_ini_type = "cosine"        ! initial temperature profile type
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



    

    !read input parameters
    call read_namelist_ex7('data/namelist/ex_7_cos.nml', nx, ny, a_adv, a_diff, total_time, max_err, Ra, T_ini_type)
    
    ! print namelist values
    print *, 'nx = ', nx
    print *, 'ny = ', ny
    print *, 'a_adv = ', a_adv
    print *, 'a_diff = ', a_diff
    print *, 'total_time = ', total_time
    print *, 'max_err = ', max_err
    print *, 'Ra = ', Ra
    print *, 'T_ini_type = ', T_ini_type

    !nx = 10
    !ny = 10
    !T_ini_type = 'empty'

    ! Init Temperature arrays
    allocate(T(nx, ny), T_new(nx, ny), w(nx, ny), psi(nx, ny), Tdx(nx, ny), u(nx, ny), v(nx, ny))
    allocate(advection(nx, ny), d2T(nx, ny))


    ! Initialize arrays
    T = 0.0
    T_new = 0.0
    w = 0.0
    psi = 0.0
    Tdx = 0.0
    u = 0.0
    v = 0.0
    d2T = 0.0
    advection = 0.0


    h = 1.0 / (ny - 1)

    nsteps = int(total_time / h)

    !call print_matrix(T)
    
    if (trim(T_ini_type) == 'rand') then
        ! Fill T with random numbers
        call RANDOM_NUMBER(T)
        T(:, 1) = 1.0
        T(:, ny) = 0.0
    else if (trim(T_ini_type) == 'cosine') then  
        call initialize_temperature_cosine(T, nx, ny)

    else if (trim(T_ini_type) == 'empty') then   
        T = 0.0
    else
        ! Handle unknown initialization type
        print *, 'Invalid initialization state:', trim(T_ini_type)
        stop
    end if

    outputfilename = 'data/ex_7/T_cos.csv'
    open(unit=10, file=trim(outputfilename), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    
    close(10)

    call write_to_csv_real8('data/ex_7/T_cos.csv', T)


    
    

    k = 0
    do while (time < total_time)
        call boundaries_T(T)

        

        ! calc first derivertiv of temperature in x direction:
        !Compute dT/dx
        call dx(T, h, Tdx)


        ! Determine 𝜔 from Ra dT/dx using the Poisson solver (Equation 2)
        f_norm = SQRT(SUM((-Ra*Tdx)**2)/(nx*ny))
        res_rms = f_norm
        do while (res_rms/f_norm > max_err)
            res_rms = Vcycle_2DPoisson(w, -Ra*Tdx, h, alpha) 
        end do


        ! Compute stream function 𝜓 from 𝜔 using the Poisson solver (Equation 3)
        f_norm = SQRT(SUM((-w)**2)/(nx*ny))
        res_rms = f_norm
        do while (res_rms/f_norm > max_err)
            res_rms = Vcycle_2DPoisson(psi, -w, h, alpha) 

        end do



        ! Compute the wind speeds 𝑢 and 𝑣 from 𝜓
        call dy(psi, h, u)
        call dx(-psi, h, v)


        ! Compute the time step from a_adv, a_diff and the maximum wind speed in the
        ! we have no kappa this time like in ex5
        dt = MIN(a_diff*h**2, a_adv*h/MAX(MAXVAL(ABS(u)), MAXVAL(ABS(v))))
        
        call FiniteDifference2D_real8(T, h, d2T)

        ! 𝑣⃑ * ∇𝑇 = advection
        call advection2D_real8(T, h, u, v, advection)
        ! calc next T
        T_new = T + dt * (d2T - advection)

        T = T_new


        ! track time and iterations
        time = time + dt
        k = k + 1
        
        ! Maxiter to protect from to long runtimes and make numerically instabel behavior debugging easier
        if (k >= 2000) then
            exit
        end if

        ! only plot a certain amount of matix to not hav to large files
        if (mod(k, 10) == 0) then
            call write_to_csv_real8('data/ex_7/T_cos.csv', T)
        end if
    end do


end program dry_convection

