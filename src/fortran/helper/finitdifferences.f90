module finitedifference
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
        

        do i = 2, n - 1
            do j = 2, m - 1
                second_derivative(i,j) = (y(i+1,j) + y(i-1,j) + y(i,j+1) - 4.0 * y(i,j) + y(i,j-1)) / h**2
            end do
        end do

    end subroutine FiniteDifference2D
end module finitedifference


