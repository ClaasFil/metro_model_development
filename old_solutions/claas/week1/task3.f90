program statistics
    implicit none
    integer :: n, i
    real :: number, sum, sum_sq, mean, std_dev

    sum = 0.0
    sum_sq = 0.0

    print *, 'Enter the number of real numbers:'
    read *, n

    do i = 1, n
        read *, number
        sum = sum + number
        sum_sq = sum_sq + number**2
    end do

    mean = sum / n
    std_dev = sqrt((sum_sq / n) - (mean**2))

    print *, 'Mean =', mean
    print *, 'St Dev =', std_dev
end program statistics
