module compute_module
    implicit none
contains
    subroutine compute(numbers, mean, std_dev)
        implicit none
        real, intent(out) :: mean, std_dev
        real, intent(in) :: numbers(:)  ! Assumed-shape array
        real :: total_sum, sum_sq
        integer :: n


        n = size(numbers)

        total_sum = sum(numbers)

        sum_sq = sum(numbers**2)

        mean = total_sum / n
        std_dev = sqrt((sum_sq / n) - (mean**2))




    end subroutine compute
end module compute_module




program statistics
    use compute_module
    implicit none
    integer :: n, i
    real :: number, mean, std_dev
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




    call compute(numbers, mean, std_dev)

    print *, 'Mean:', mean
    print *, 'Standard deviation:', std_dev

    deallocate(numbers)
end program statistics










