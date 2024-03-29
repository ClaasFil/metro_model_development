program statistics
    implicit none
    integer ::  i, j, n
    real ::  h
    real, allocatable :: arrIn(:), matrix(:,:), arrOut(:), analytical_solution(:)



    ! Ask for the number of grid points and grid spacing
    print *, "Enter the number of grid points (n):"
    read(*, *) n
    print *, "Enter the grid spacing (h):"
    read(*, *) h



    ! Allocate arrays
    allocate(arrIn(n), analytical_solution(n))




    ! Initialize y array with some example function (e.g., sin(x), x**2)
    do i = 1, n
        arrIn(i) = (real(i)*h)**4! Example function: x**4
        analytical_solution(i) = 12 * (real(i)*h) ** 2 ! (real(i) * h)! Second derivative of example function: x**4
        !arrIn(i) = sin(real(i) * h) ! Example function: sin(x)
        !analytical_solution(i) = -sin(real(i) * h) ! Second derivative of example function: sin(x)
    end do




    ! --------- creat mtx for the 2nd derivative ---------
    ! calculate the 2nd derivative via matrix
    ! creating the matrix if quit expencive but one time cost
    ! calculating the 2nd derivative is cheap and usable many times
    
    allocate(matrix(n, n))

    ! Initialize the matrix with 0s
    matrix = 0.0

    ! Set -2 on the diagonal and 1 above and below the diagonal
    do i = 2, n -1 
        matrix(i,i) = -2.0    ! Set -2 on the diagonal
        if (i < n) then
            matrix(i, i +1) = 1.0 ! Set 1 above the diagonal
        end if
        if (i > 1) then
            matrix(i, i - 1) = 1.0 ! Set 1 above the diagonal
        end if
    end do

    !Print the matrix
    !print *, 'Matrix:'
    !do i = 1, n
    !    print *, (matrix(i,j), j = 1, n)
   ! end do



    !------ calc 2nd derivative of arrIn and store in arrOut

    allocate(arrOut(n))
    arrOut = matmul(matrix, arrIn)/h**2


    print *, 'numerical 2nd derivative:'
    print *, arrOut

    ! Output the analytical second derivative
    print *, "Second derivative analytically calculated:"
    print *, analytical_solution

    ! print error:
    print *, "Error:"
    print *, arrOut - analytical_solution



    !save to file
    open(1, file='task3_res.txt', status='replace')
    write(1, *) 'numerical 2nd derivative:'
    write(1, *) arrOut
    write(1, *) "Second derivative analytically calculated:"
    write(1, *) analytical_solution
    write(1, *) "Error:"
    write(1, *) arrOut - analytical_solution


    ! Deallocate the matrix
    deallocate(matrix)
    
    end program statistics
