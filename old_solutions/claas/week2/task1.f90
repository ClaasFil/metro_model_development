program statistics
    implicit none
    integer :: n, i
    real :: number, total_sum, mean, sum_sq, std_dev
    real, allocatable :: numbers(:)

    

    print *, 'Enter the number of real numbers:'
    read *, n

    do 
        if (n > 0) then
            exit
        end if
        print *, 'Invalid input enter positiv number'
        read *, n
    end do


    allocate(numbers(n))

    print *, 'Enter the real numbers:'
    do i = 1, n
        read *, number
        numbers(i) = number
    end do

    
    total_sum = sum(numbers)

    sum_sq = sum(numbers**2)

    mean = total_sum / n
    std_dev = sqrt((sum_sq / n) - (mean**2))

    print *, 'Mean =', mean
    print *, 'St Dev =', std_dev
end program statistics
