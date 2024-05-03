program HeatDiffusion2d
    use finitedifference
    use namelist_utilities
    use csv_writer
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx, ny, i, j, k, nsteps, io_error
    real :: L, total_time, h, dt, time, a, kappa
    real, allocatable :: T_init(:,:), T(:,:), T_new(:,:), second_derivative(:,:)
    character(len=30) :: outputfilename

    ! Initialise Parameters
    L = 1.0           ! Length of the domain
    nx = 64
    ny = 64
    kappa = 1.0
    total_time = 0.1
    a = 0.3
    outputfilename='data/heatdif2d_out_random_a03'


    call read_namelist('data/namelist/namelist.nml', nx, ny, kappa, a, outputfilename)

    ! Calculate grid spacing and time step
    h = L / real(nx - 1)
    dt = a * h**2 / kappa
    nsteps = int(total_time / dt)

    ! Allocate memory for temperature arrays
    allocate(T_init(nx, ny), T(nx, ny), T_new(nx, ny), second_derivative(nx, ny))

    ! Initialize temperature array
    T_init = 0.0
    T_init(nx/2, ny/2) = 10.0  ! Spike in the center
    call RANDOM_NUMBER(T_init)
    T = T_init

    outputfilename = 'data/heatdif2d_out.csv'
    open(unit=10, file=trim(outputfilename), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif

    write(10, '(A, I0, A, I0)') 'nx,', nx, ',ny,', ny
    close(10)

    ! Integration loop
    do k = 1, int(nsteps)
        ! Compute second derivative
        call FiniteDifference2D(T, h, second_derivative)
        
        ! Update temperature using forward Euler method
        T_new = T + dt * kappa * second_derivative
        
        ! Apply boundary conditions
        T_new(1, 1) = 0.0
        T_new(nx, ny) = 0.0

        ! Update temperature array for next time step
        T = T_new

        call write_to_csv(outputfilename, T)

    end do

    close(10)

    ! Deallocate memory
    deallocate(T_init, T, T_new, second_derivative)

    
end program HeatDiffusion2d











