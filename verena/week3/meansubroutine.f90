program MeanSubroutine

    implicit none
    integer :: n, i
    real, allocatable :: numbers(:)
    real :: mean, std_deviation

    ! Main program
    print *, "How many numbers do you want to enter?"
    read(*, *) n

    ! Check if number is negative
    do while (n <= 0)
        print *, "Error: Please enter a positive number!"
        read(*, *) n
    end do

    allocate(numbers(n))

    print *, "Enter the numbers one by one: "
    read(*, *) numbers

    ! Call subroutine to compute mean and standard deviation
    call compute_mean_and_std(numbers, mean, std_deviation)

    print *, "Mean: ", mean
    print *, "Standard Deviation: ", std_deviation

    ! Subroutine
    contains
        subroutine compute_mean_and_std(array, mean, std_deviation)
            real, intent(in) :: array(:)
            real, intent(out) :: mean, std_deviation
            real :: sum_numbers, sum_of_squares
            integer :: n

            n = size(array)
            sum_numbers = sum(array)
            sum_of_squares = sum(array**2)
            mean = sum_numbers / real(n)
            std_deviation = sqrt(sum_of_squares / real(n) - mean**2)
        end subroutine compute_mean_and_std

end program MeanSubroutine