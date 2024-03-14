program factorial
    implicit none
    integer :: n, i
    integer, dimension(:), allocatable :: fact

    print *, 'Enter a positive integer:'
    read *, n

    if (n < 1) then
        print *, 'Error: The number must be positive.'
        stop
    end if

    allocate(fact(n))
    fact(1) = 1
    do i = 2, n
        fact(i) = i * fact(i-1)
    end do

    print *, 'The factorial of ', n, ' is ', fact(n)

    deallocate(fact)
end program factorial