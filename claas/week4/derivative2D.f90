
module matrix_utilities
    implicit none
contains
    subroutine print_matrix(matrix)
        real, dimension(:, :), intent(in) :: matrix
        integer :: i, j, n, m

        n = size(matrix, 1)
        m = size(matrix, 2)

        ! Print the matrix in a formatted manner
        print *, 'Matrix:'
        do i = 1, n
            ! Print each row with aligned columns
            write(*, '(100F8.2)') (matrix(i, j), j = 1, m)
        end do
    end subroutine print_matrix
end module matrix_utilities



module derivative2D
    implicit none
    
contains
    subroutine deriv2(y, h, second_derivative)
        use matrix_utilities
        real, allocatable, intent(in) :: y(:,:)
        real, intent(in) :: h
        real, allocatable, intent(out) :: second_derivative(:,:)
        integer :: n, m, i, j
        
        n = size(y, 1)
        m = size(y, 2)
        
     
        
        

        allocate(second_derivative(n,m))

        second_derivative = 0.0

        do i = 2, n - 1
            do j = 2, m - 1
  
                second_derivative(i,j) = (y(i+1,j) + y(i-1,j) + y(i,j+1) - 4.0 * y(i,j) + y(i,j-1)) / h**2
            end do
        end do

  
      
      


        
    end subroutine deriv2
end module derivative2D














