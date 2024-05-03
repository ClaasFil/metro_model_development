program SieveOfEratosthenes
    implicit none
    integer :: n, i, j
    logical, allocatable :: is_prime(:)

    ! Ask for maximum number
    print *, "Enter the maximum number to find primes up to:"
    read(*, *) n

    ! Allocate an array to keep track of whether each number is prime
    allocate(is_prime(n))
    is_prime = .true.

    is_prime(1) = .false. ! 1 is not prime
    is_prime(2) = .true.  ! 2 is prime

    ! Apply Sieve of Eratosthenes algorithm
    do i = 2, n
        if (is_prime(i)) then
            ! Mark multiples of i as not prime
            do j = 2*i, n, i
                is_prime(j) = .false.
            end do
        end if
    end do

    ! Output prime numbers
    print *, "Prime numbers up to", n, ":"
    do i = 2, n
        if (is_prime(i)) then
            print *, i
        end if
    end do

    deallocate(is_prime)

end program SieveOfEratosthenes
