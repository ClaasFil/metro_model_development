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
    integer :: nx=257, ny=65                    ! number of grid points in x and y
    integer :: io_error                      ! error code for file IO
    integer :: k                            ! iteration counter
    integer :: nsteps                        ! number of time steps
    integer :: i, j                          ! loop variables
    real :: a_adv = 0.4                 ! advection coefficients 
    real :: a_diff = 0.23                   ! diffusion coefficients
    real :: total_time = 0.1                      ! total simulation time
    real :: max_err = 1.E-3                          ! maximum error
    real :: Ra = 1.E5                               ! Rayleigh number
    real :: Pr = 1                             ! Prandtl number
    real :: gamma = 0.5                 ! Dry adiabatic lapse rate
    real :: lambda = 0.2                    ! Latent heating
    real :: tau = 2.0                      ! Time scale for condensation
    real :: alpha = 1.0                       !  Exponential scaling factor for saturation specific humidity
    real :: Sm = 1.33                       ! Ratio of diffusivities
    real(8) :: h                             ! grid spacing
    real(8) :: dt                            ! time step
    real(8) :: relax = 1.                         ! relaxation parameter
    REAL(8) :: res_rms                          ! root mean square residue 
    REAL(8) :: time = 0.0                     ! time. since h is change wee neet to track time
    REAL(8) :: f_norm                          ! norm used for convergence check of poisson solver
    character(len=100) :: T_ini_type = "cosine"        ! initial temperature profile type
    character(len=100) :: outputfilename_T        ! output file name
    character(len=100) :: outputfilename_q        ! output file name
    character(len=100) :: outputfilename_C        ! output file name
    character(len=100) :: outputfilename_psi        ! output file name
    real(8), allocatable:: T(:,:)               ! temperature profile
    real(8), allocatable:: T_new(:,:)           ! new temperature profile
    real(8), allocatable:: w(:,:)               ! vorticity
    real(8), allocatable:: w_new(:,:)           ! new vorticity
    real(8), allocatable:: psi(:,:)             ! stream function
    real(8), allocatable:: Tdx(:,:)             ! temperature derivative in x
    real(8), allocatable:: u(:,:)               ! x velocity
    real(8), allocatable:: v(:,:)               ! y velocity
    real(8), allocatable:: d2T(:,:)             ! second derivative(x,y) of temperature
    real(8), allocatable:: d2w(:,:)             ! second derivative(x,y) of vorticity
    real(8), allocatable:: d2q(:,:)             ! second derivative(x,y) of specific humidity
    real(8), allocatable:: vdT(:,:)             ! advection
    real(8), allocatable:: vdw(:,:)             ! advection
    real(8), allocatable:: vdq(:,:)             ! advection
    real(8), allocatable:: RadTdx(:,:)          ! Ra * dT/dx
    real(8), allocatable :: q(:,:)                        ! specific humidity
    real(8), allocatable :: q_new(:,:)                    ! specific humidity
    real(8), allocatable :: qs(:,:)                    ! saturation specific humidity
    real(8), allocatable :: C(:,:)                    ! condensation
    real(8), allocatable :: q_qs(:,:)                    ! heaviside function
    

    !read input parameters
    call read_namelist_ex9('data/namelist/ex_9_cos.nml', nx, ny, a_adv, a_diff, total_time, max_err, Ra, &
                            T_ini_type, Pr, alpha, gamma, lambda, tau)
    
    ! print namelist values
    print *, 'nx = ', nx
    print *, 'ny = ', ny
    print *, 'a_adv = ', a_adv
    print *, 'a_diff = ', a_diff
    print *, 'total_time = ', total_time
    print *, 'max_err = ', max_err
    print *, 'Ra = ', Ra
    print *, 'T_ini_type = ', T_ini_type
    print *, 'Pr = ', Pr
    print *, 'alpha = ', alpha
    print *, 'gamma = ', gamma
    print *, 'lambda = ', lambda
    print *, 'tau = ', tau

    !nx = 10
    !ny = 10
    !T_ini_type = 'empty'

    ! Init Temperature arrays
    allocate(T(nx, ny), T_new(nx, ny), w(nx, ny), w_new(nx, ny), psi(nx, ny), Tdx(nx, ny), u(nx, ny), v(nx, ny))
    allocate(vdT(nx, ny), vdw(nx,ny), vdq(nx, ny), d2T(nx, ny), d2w(nx, ny), d2q(nx, ny), RadTdx(nx, ny))
    allocate(q(nx, ny), qs(nx, ny), q_new(nx, ny), C(nx, ny), q_qs(nx, ny))


    ! Initialize arrays
    T = 0.0
    T_new = 0.0
    w = 0.0
    w_new = 0.0
    psi = 0.0
    Tdx = 0.0
    u = 0.0
    v = 0.0
    d2T = 0.0
    d2w = 0.0

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

    qs = exp(alpha * T)
    q = 0.8 * qs

    outputfilename_T = 'data/ex_9/T_cos.csv'
    outputfilename_q = 'data/ex_9/q_cos.csv'
    outputfilename_C = 'data/ex_9/C_cos.csv'
    outputfilename_psi = 'data/ex_9/psi_cos.csv'
    open(unit=10, file=trim(outputfilename_T), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    open(unit=10, file=trim(outputfilename_q), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    open(unit=10, file=trim(outputfilename_C), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)

    open(unit=10, file=trim(outputfilename_psi), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    close(10)


    call write_to_csv_real8(outputfilename_T, T)
    call write_to_csv_real8(outputfilename_q, q)
    call write_to_csv_real8(outputfilename_C, C)
    call write_to_csv_real8(outputfilename_psi, psi)


    k = 0
    do while (time < total_time)

        call boundaries_T(T)

        ! Compute stream function 𝜓 from 𝜔 using the Poisson solver (Equation 2)
        f_norm = SQRT(SUM((-w)**2)/(nx*ny))
        res_rms = f_norm
        do while (res_rms/f_norm > max_err)
            res_rms = Vcycle_2DPoisson(psi, -w, h, relax) 
        end do

        ! Compute the wind speeds 𝑢 and 𝑣 from 𝜓
        call dy(psi, h, u)
        call dx(-psi, h, v)

        ! Compute the time step from a_adv, a_diff and the maximum wind speed in the domain
        dt = MIN(a_diff*h**2/max(1.,Pr), a_adv*h/MAX(MAXVAL(ABS(u)), MAXVAL(ABS(v))))  

        ! Calc first derivertiv of temperature in x direction:
        ! Compute Ra * dT/dx
        call dx(T, h, Tdx)
        RadTdx = Ra * Tdx

        ! Compute the second derivative of temperature
        call FiniteDifference2D_real8(T, h, d2T)
        ! 𝑣⃑ * ∇𝑇 = advection
        call advection2D_real8(T, h, u, v, vdT)

        ! Compute the second derivative of vorticity 𝜔
        call FiniteDifference2D_real8(w, h, d2w)
        ! 𝑣⃑ * ∇𝜔 = advection
        call advection2D_real8(w, h, u, v, vdw)

        ! Compute the second derivative of specific humdity q
        call FiniteDifference2D_real8(q, h, d2q)
        ! 𝑣⃑ * ∇q = advection
        call advection2D_real8(q, h, u, v, vdq)

        ! calc next T
        T_new = T + dt * (d2T - vdT - gamma * v) 
        T = T_new

        ! Calc next q
        q_new = q + dt * (Sm*d2q - vdq) 
        q = q_new

        ! calc next 𝜔
        w_new = w + dt * (Pr * d2w + Pr * RadTdx - vdw)
        w = w_new
        call boundaries_zero(w)

        ! track time and iterations
        time = time + dt
        k = k + 1

        !Diagnose condensation
        q_qs = q - qs
        !call print_matrix(q_qs)
        C = (q_qs) / tau * heaviside(q_qs)
        !print *, 'Maximu of C = ', MAXVAL(C)
        print * , "C>0 ", COUNT(C > 0.0)
        !print *, 'Minimu of C = ', MINVAL(C)

        ! Update T and q
        T_new = T + lambda * C
        T = T_new
        q_new = q - C
        q = q_new

        
        ! Maxiter to protect from to long runtimes and make numerically instable behavior debugging easier
        if (k >= 2500) then
            exit
        end if

        ! only plot a certain amount of matrix to not have to large files
        !if (mod(k, 10) == 0) then
            !call write_to_csv_real8(outputfilename, T)
        !end if
    end do
    call write_to_csv_real8(outputfilename_T, T)
    call write_to_csv_real8(outputfilename_q, q)
    call write_to_csv_real8(outputfilename_C, C)
    call write_to_csv_real8(outputfilename_psi, psi)

    contains

        FUNCTION heaviside(x) RESULT(y)
            REAL(8), INTENT(IN), allocatable :: x(:,:)
            REAL(8), allocatable :: y(:,:)

            nx = size(x, 1)
            ny = size(x, 2)

            allocate(y(nx, ny))

            do CONCURRENT(i = 1:nx, j = 1:ny)
                if (x(i, j) > 0.0) then
                    y(i, j) = 1.0
                    !print *, 'x(i, j) = ', x(i, j)
                else
                    y(i, j) = 0.0
                end if
            end do
            
        END FUNCTION heaviside

        FUNCTION qsat(T) RESULT(q)
            REAL(8), INTENT(IN), allocatable:: T(:,:)
            REAL(8), allocatable :: q(:,:)

            nx = size(T, 1)
            ny = size(T, 2)
            allocate(q(nx, ny))

            q = EXP(alpha * T)
        END FUNCTION qsat

        


end program dry_convection

