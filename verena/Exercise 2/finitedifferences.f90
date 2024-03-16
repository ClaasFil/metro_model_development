program FiniteDifferences
    implicit none
    integer :: n, i
    real, allocatable :: y(:), second_derivative(:)
    real :: h

    ! Ask for the number of grid points and grid spacing
    print *, "Enter the number of grid points (n):"
    read(*, *) n
    print *, "Enter the grid spacing (h):"
    read(*, *) h

    ! Allocate arrays
    allocate(y(n), second_derivative(n), analytical_solution(n))

    ! Initialize y array with some example function (e.g., sin(x))
    do i = 1, n
        y(i) = sin(real(i) * h) ! Example function: sin(x)
        analytical_solution(i) = -sin(real(i) * h) ! Second derivative of example function: sin(x)
    end do

    ! Compute second derivative using finite difference method
    second_derivative(1) = 0.0
    do i = 2, n - 1
        second_derivative(i) = (y(i+1) - 2.0 * y(i) + y(i-1)) / h**2
    end do
    second_derivative(n) = 0.0

    ! Output the second derivative
    print *, "Second derivative using finite difference method:"
    do i = 1, n
        print *, second_derivative(i)
    end do

    ! Output the analytical second derivative
    print *, "Second derivative analytically calculated:"
    do i = 1, n
        print *, analytical_solution(i)
    end do

    ! Output the second derivative
    print *, "Error:"
    do i = 1, n
        print *, abs(second_derivative(i)-analytical_solution(i))
    end do

    deallocate(y, second_derivative, analytical_solution)

end program FiniteDifferences