module finitedifference
    use matrix_utilities
    implicit none
    
contains
    subroutine FiniteDifference2D(y, h, second_derivative)
        real, intent(in) :: y(:,:)
        real, intent(in) :: h
        real, intent(out) :: second_derivative(:,:)
        integer :: n, m, i, j
        
        n = size(y, 1)
        m = size(y, 2)
        

        second_derivative = 0.0
        

        do j = 2, m - 1
            do i = 2, n - 1
                second_derivative(i,j) = (y(i+1,j) + y(i-1,j) + y(i,j+1) - 4.0 * y(i,j) + y(i,j-1)) / h**2
            end do
        end do

        ! Apply boundary conditions
        ! Top and bottom rows
        do j = 1, m
            second_derivative(1,j) = (y(2,j) - y(1,j)) / h**2  ! Top row
            second_derivative(n,j) = (y(n-1,j) - y(n,j)) / h**2  ! Bottom row
        end do

        ! Left and right columns
        do i = 2, n - 1
            second_derivative(i,1) = (y(i,2) - y(i,1)) / h**2  ! Left column
            second_derivative(i,m) = (y(i,m-1) - y(i,m)) / h**2  ! Right column
        end do

    end subroutine FiniteDifference2D

    subroutine FiniteDifference2D_real8(y, h, second_derivative)
        real(8), intent(in) :: y(:,:)
        real(8), intent(in) :: h
        real(8), intent(out) :: second_derivative(:,:)
        integer :: n, m, i, j
        
        !call print_matrix(y)
        !print *, 'h' , h

        n = size(y, 1)
        m = size(y, 2)

        second_derivative = 0.0
        

        do j = 2, m - 1
            do i = 2, n - 1
                second_derivative(i,j) = 1./h**2 * (y(i,j+1) + y(i,j-1) + y(i+1,j) + y(i-1,j) - 4*y(i,j))
            end do
        end do

        ! Apply boundary conditions
        ! Top and bottom rows
        !do j = 1, m
        !    second_derivative(1,j) = (y(2,j) - y(1,j)) / h**2  ! Top row
        !    second_derivative(n,j) = (y(n-1,j) - y(n,j)) / h**2  ! Bottom row
        !end do

        ! Left and right columns
        !do i = 1, n
        !    second_derivative(i,1) = (y(i,2) - y(i,1)) / h**2  ! Left column
        !    second_derivative(i,m) = (y(i,m-1) - y(i,m)) / h**2  ! Right column
        !end do

    end subroutine FiniteDifference2D_real8


    subroutine dx(y, h, first_derivative)
        ! Compute the first derivative of a 2D matrix in the x direction
        real(8), intent(in) :: y(:,:)   ! Input matrix
        real(8), intent(in) :: h        ! Grid spacing
        real(8), intent(out) :: first_derivative(:,:)  ! Output matrix for the first derivative
        integer :: n, m, i, j
        
        n = size(y, 1)
        m = size(y, 2)
        
        first_derivative = 0.0  
        
        do j = 1, m
            do i = 2, n - 1
                first_derivative(i, j) = (y(i + 1, j) - y(i - 1, j)) / (2 * h)  ! Centered difference
            end do
        end do

    end subroutine dx


    subroutine dy(y, h, first_derivative)
        ! Compute the first derivative of a 2D matrix in the y direction
        real(8), intent(in) :: y(:,:)   ! Input matrix
        real(8), intent(in) :: h        ! Grid spacing
        real(8), intent(out) :: first_derivative(:,:)  ! Output matrix for the first derivative
        integer :: n, m, i, j
        
        n = size(y, 1)
        m = size(y, 2)
        
        first_derivative = 0.0  
        
        do j = 2, m - 1
            do i = 1, n
                first_derivative(i, j) = (y(i, j + 1) - y(i, j - 1)) / (2 * h)  ! Centered difference
            end do
        end do

    end subroutine dy





    subroutine advection2D(T, h, u, v, advection)
        real, intent(in) :: T(:,:)
        real, intent(in) :: h
        real, intent(in) :: u(:,:)
        real, intent(in) :: v(:,:)
        real, intent(out) :: advection(:,:)
        integer :: n, m, i, j
        
        n = size(T, 1)
        m = size(T, 2)
        

        advection = 0.0
        

        do j = 2, m - 1
            do i = 2, n - 1
                ! Calculate advection term using upwind differencing
                ! For u component
                if (u(i,j) >= 0.0) then
                    advection(i,j) = advection(i,j) + u(i,j) * (T(i,j) - T(i-1,j)) / h
                else
                    advection(i,j) = advection(i,j) + u(i,j) * (T(i+1,j) - T(i,j)) / h
                end if
                
                ! For v component
                if (v(i,j) >= 0.0) then
                    advection(i,j) = advection(i,j) + v(i,j) * (T(i,j) - T(i,j-1)) / h
                else
                    advection(i,j) = advection(i,j) + v(i,j) * (T(i,j+1) - T(i,j)) / h
                end if
            end do
        end do

        !call print_matrix(advection)

    end subroutine advection2D


    subroutine advection2D_real8(T, h, u, v, advection)
        real(8), intent(in) :: T(:,:)
        real(8), intent(in) :: h
        real(8), intent(in) :: u(:,:)
        real(8), intent(in) :: v(:,:)
        real(8), intent(out) :: advection(:,:)
        integer :: n, m, i, j
        
        n = size(T, 1)
        m = size(T, 2)
        

        advection = 0.0
        

        do j = 2, m - 1
            do i = 2, n - 1
                ! Calculate advection term using upwind differencing
                ! For u component
                if (u(i,j) >= 0.0) then
                    advection(i,j) = advection(i,j) + u(i,j) * (T(i,j) - T(i-1,j)) / h
                else
                    advection(i,j) = advection(i,j) + u(i,j) * (T(i+1,j) - T(i,j)) / h
                end if
                
                ! For v component
                if (v(i,j) >= 0.0) then
                    advection(i,j) = advection(i,j) + v(i,j) * (T(i,j) - T(i,j-1)) / h
                else
                    advection(i,j) = advection(i,j) + v(i,j) * (T(i,j+1) - T(i,j)) / h
                end if
            end do
        end do

        !call print_matrix(advection)

    end subroutine advection2D_real8
end module finitedifference


