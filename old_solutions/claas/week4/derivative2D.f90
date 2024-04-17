
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

module temperature_io
    implicit none
contains
    subroutine write_temperature_to_csv(outfile, T)
        implicit none
        real, dimension(:,:), intent(in) :: T
        character(len=*), intent(in)  :: outfile
        integer :: io, i, j

        ! Open the file for appending; create a new one if it doesn't exist
        open(unit=20, file=trim(outfile), status='unknown', action='write', iostat=io, position='append')
        if (io /= 0) then
            print *, "Failed to open file:", trim(outfile)
            return
        end if


        ! Write the matrix to the file
        do i = 1, size(T, 1)
            do j = 1, size(T, 2)
                if (j == size(T, 2)) then
                    write(20, '(F6.2)') T(i, j)  ! Last element in the row
                else
                    write(20, '(F6.2, A)', advance='no') T(i, j), ','  ! Elements with comma
                end if
            end do
            write(20, *)  ! Newline for the next row
        end do

        ! Close the file
        close(20)
    end subroutine write_temperature_to_csv
end module temperature_io

module derivative2D
    implicit none
    
contains
    subroutine deriv2(y, h, second_derivative)
        use matrix_utilities
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

  
      
      


        
    end subroutine deriv2
end module derivative2D





module norm_calculations
    implicit none
contains
    subroutine calc_linfty_norm(matrix, norm_linfty)
        ! Declare the subroutine arguments
        real, intent(in) :: matrix(:,:)
        real, intent(out) :: norm_linfty
        integer :: n, m, i, j
        real :: row_sum

        ! Get the dimensions of the matrix
        n = size(matrix, 1)
        m = size(matrix, 2)

        ! Initialize the L-infinity norm
        norm_linfty = 0.0

        ! Compute L-infinity norm (maximum row sum)
        do i = 1, n
            row_sum = 0.0
            do j = 1, m
                
                row_sum = row_sum + abs(matrix(i, j))
            end do
            norm_linfty = max(norm_linfty, row_sum)
        end do
    end subroutine calc_linfty_norm
end module norm_calculations








