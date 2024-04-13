program HeatDiffusion2d
    use FiniteDifference2d
    implicit none
    integer :: nx, ny, i, j, k, nsteps
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
    outputfilename='heatdif2d_out_random_a03.txt'

    ! Read in Parameters
    NAMELIST /inputs/ nx, ny, kappa, total_time, a, outputfilename
    open(1, file="namelist", status="old")
    read(1, NML=inputs)
    close(1)
    write(*, NML=inputs)

    ! Read in Parameters
    ! open(1, file="namelist.txt", status="old")
    ! read(1, *) nx, ny, kappa, total_time, a, outputfilename
    ! close(1)
    ! write(*, *) nx, ny, kappa, total_time, a, outputfilename

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

    ! Write out temperature
    open(unit=10, file=outputfilename, status='unknown')
    write(10, *) T

    ! Integration loop
    do k = 1, int(nsteps)
        ! Compute second derivative
        call finitedifferences(T, h, second_derivative)
        
        ! Update temperature using forward Euler method
        do i = 2, nx - 1
            do j = 2, ny - 1
                T_new(i, j) = T(i, j) + dt * kappa * second_derivative(i, j)
            end do
        end do
        
        ! Apply boundary conditions
        T_new(1, 1) = 0.0
        T_new(nx, ny) = 0.0

        ! Update temperature array for next time step
        T = T_new
        if (mod(k, 10) == 0) then
            write(10, *) T
        end if
    end do

    close(10)

    ! Deallocate memory
    deallocate(T_init, T, T_new, second_derivative)

end program HeatDiffusion2d