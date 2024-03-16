program statistics
    implicit none
    integer :: h , i , j
    real :: number 
    real, allocatable :: arrIn(:), matrix(:,:), arrOut(:)




    print *, 'Enter the number of real numbers:'
    read *, h

    do 
        if (h > 0) then
            exit
        end if
        print *, 'Invalid input enter positiv number'
        read *, h
    end do


    allocate(arrIn(h))

    print *, 'Enter the real numbers:'
    do i = 1, h
        read *, number
        arrIn(i) = number
    end do


    ! --------- creat mtx for the 2nd derivative ---------
    ! calculate the 2nd derivative via matrix
    !creating the matrix if quit expencive but one time cost
    ! calculating the 2nd derivative is cheap and usable many times
    
    allocate(matrix(h, h))

    ! Initialize the matrix with 0s
    matrix = 0.0

    ! Set -2 on the diagonal and 1 above and below the diagonal
    do i = 2, h -1 
        matrix(i,i) = -2.0    ! Set -2 on the diagonal
        if (i < h) then
            matrix(i, i +1) = 1.0 ! Set 1 above the diagonal
        end if
        if (i > 1) then
            matrix(i, i - 1) = 1.0 ! Set 1 above the diagonal
        end if
    end do

    ! Print the matrix
    !print *, 'Matrix:'
    !do i = 1, h
    !    print *, (matrix(i,j), j = 1, h)
    !end do



    !------ calc 2nd derivative of arrIn and store in arrOut

    allocate(arrOut(h))
    arrOut = matmul(matrix, arrIn)


    print *, '2nd derivative:'
    print *, arrOut


    ! Deallocate the matrix
    deallocate(matrix)
    
    end program statistics
