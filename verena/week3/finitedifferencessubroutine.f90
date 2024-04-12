module finitedifferencessubroutine
    implicit none
    
contains
    subroutine finitedifferences(y, h, second_derivative)
        real, allocatable, intent(in) :: y(:)
        real, intent(in) :: h
        real, allocatable, intent(out) :: second_derivative(:)
        integer :: n, i

        n = size(y)

        allocate(second_derivative(n))

        second_derivative(1) = 0.0
        do i = 2, n - 1
            second_derivative(i) = (y(i+1) - 2.0 * y(i) + y(i-1)) / h**2
        end do
        second_derivative(n) = 0.0

        return
    end subroutine finitedifferences
end module finitedifferencessubroutine
