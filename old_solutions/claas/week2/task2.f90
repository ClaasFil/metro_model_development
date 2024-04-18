program statistics
    implicit none
    integer :: max_number, index, each_Number
    logical, allocatable :: numbers(:)


    

    print *, 'Enter the number of real numbers:'
    read *, max_number

    do 
        if (max_number > 0) then
            exit
        end if
        print *, 'Invalid input enter positiv number'
        read *, max_number
    end do


    allocate(numbers(max_number))

    !print *, numbers

    do index = 2, max_number
        if (numbers(index) .eqv. .false.) then
            !print *, 'Index: ', index    
            !print *, numbers
            do each_Number = index + 1 , max_number
                
                if (mod(each_Number, index) == 0) then
                    numbers(each_Number) = .true.
                    !print *, each_Number
                end if
            end do
        end if
        

    end do

    !print *, numbers

! return al indet which are false
    do index = 2, max_number
        if (numbers(index) .eqv. .false.) then
            print *, index
        end if
    end do

    deallocate(numbers)


end program statistics
