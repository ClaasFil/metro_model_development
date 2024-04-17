program HeatDiffusion
    use FiniteDifferencesSubroutine
    implicit none
    integer :: nx, i, j
    real :: L, total_time, h, dt, time
    real, allocatable :: T_old(:), T(:), T_new(:), second_derivative(:)

    ! Parameters
    L = 1.0           ! Length of the domain
    nx = 100          ! Number of grid points
    total_time = 0.1  ! Total integration time

    ! Calculate grid spacing and time step
    h = L / real(nx - 1)
    dt = 0.4 * h**2

    ! Allocate memory for temperature arrays
    allocate(T_old(nx), T_new(nx), second_derivative(nx))

    ! Initialize temperature array
    T_old = 0.0
    T_old(nx/2) = 10.0  ! Spike in the center

    T = T_old

    ! Integration loop
    do j = 1, int(total_time / dt)
        ! Compute second derivative
        call finitedifferences(T, h, second_derivative)
        
        ! Update temperature using forward Euler method
        do i = 2, nx - 1
            T_new(i) = T(i) + dt * second_derivative(i)
        end do
        
        ! Apply boundary conditions
        T_new(1) = 0.0
        T_new(nx) = 0.0

        ! Update temperature array for next time step
        T = T_new
    end do

    ! Write out temperature at the first and last time step
    open(unit=10, file='heatdiffusion_output.txt', status='unknown')
    write(10, *) T_old
    write(10, *) T
    close(10)

    ! Deallocate memory
    deallocate(T_old, T, T_new, second_derivative)

end program HeatDiffusion
