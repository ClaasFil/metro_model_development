module finitedifference
    implicit none
    
contains
    subroutine FiniteDifference2D(T, h, second_derivative)
        real, intent(in) :: T(:,:)
        real, intent(in) :: h
        real, intent(out) :: second_derivative(size(T,1), size(T,2))
        integer :: i, j

        ! Calculate second derivative in the x-direction
        do j = 1, size(T, 2)
            second_derivative(1, j) = 0.0
            do i = 2, size(T, 1) - 1
                second_derivative(i, j) = (T(i+1, j) - 2.0 * T(i, j) + T(i-1, j)) / h**2
            end do
            second_derivative(size(T, 1), j) = 0.0
        end do

        ! Calculate second derivative in the y-direction
        do i = 1, size(T, 1)
            second_derivative(i, 1) = 0.0
            do j = 2, size(T, 2) - 1
                second_derivative(i, j) = second_derivative(i, j) + (T(i, j+1) - 2.0 * T(i, j) + T(i, j-1)) / h**2
            end do
            second_derivative(i, size(T, 2)) = 0.0
        end do

    end subroutine FiniteDifference2D
end module finitedifference