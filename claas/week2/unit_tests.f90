program test_statistics
    implicit none
    integer, parameter :: h = 50  ! Example size for the input array
    real, dimension(h) :: test_array, derivative_output
    real, dimension(h,h) :: derivative_matrix
    integer :: i

    ! Define a test input array (for simplicity, use a linear function whose second derivative is known)
    test_array = [(i, i=1,h)]

    ! Expected output for the second derivative of a linear function should be zeros
    ! (since the second derivative of a linear function is 0)

    ! Call the function/subroutine from your 'statistics' program that creates the matrix
    ! Assuming it's named 'create_derivative_matrix' and it returns the matrix
    ! (You'll need to modify your 'statistics' program to include this function)
    derivative_matrix = create_derivative_matrix(h)

    ! Call the function/subroutine that calculates the second derivative
    ! Assuming it's named 'calculate_second_derivative' and it returns the derivative array
    derivative_output = calculate_second_derivative(derivative_matrix, test_array)

    ! Check the output against the expected result (zeros)
    do i = 1, h
        if (abs(derivative_output(i)) > 1.0E-6) then
            print *, 'Test failed for element', i, ': Expected 0, got', derivative_output(i)
        else
            print *, 'Test passed for element', i
        end if
    end do

end program test_statistics
