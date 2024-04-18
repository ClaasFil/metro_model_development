!gfortran -o task2 derivative.f90 task2.f90 && ./task2

program task2
    use derivative
    implicit none
    integer ::  i, j, n
    real ::  h
    real, allocatable :: arrIn_sin(:), arrIn_poly(:), matrix(:,:), &
                        arrOut_sin(:), arrOut_poly(:), analytical_solution_sin(:),&
                         analytical_solution_poly(:)



    ! Ask for the number of grid points and grid spacing
    !print *, "Enter the number of grid points (n):"
    !read(*, *) n
    !print *, "Enter the grid spacing (h):"
    !read(*, *) h

    n = 10
    h = 0.1


    ! Allocate arrays
    allocate(arrIn_sin(n), arrIn_poly(n), analytical_solution_poly(n),&
             analytical_solution_sin(n), matrix(n,n), arrOut_sin(n), &
             arrOut_poly(n))



    ! Initialize y array with some example function (e.g., sin(x), x**2)
    do i = 1, n
        arrIn_sin(i) = (real(i)*h)**4! Example function: x**4
        analytical_solution_sin(i) = 12 * (real(i)*h) ** 2 ! (real(i) * h)! Second derivative of example function: x**4
    end do



    ! Initialize y array with some example function (e.g., sin(x), x**2)
    do i = 1, n
        arrIn_poly(i) = sin(real(i) * h) ! Example function: sin(x)
        analytical_solution_poly(i) = -sin(real(i) * h) ! Second derivative of example function: sin(x)
    end do

    call deriv2(arrIn_poly, arrOut_poly, h)
    call deriv2(arrIn_sin, arrOut_sin, h)





    print *, "Analytical solution for sin(x):"
    print *, analytical_solution_sin
    print *, "numerical solution for sin(x):"
    print *, arrOut_sin
    print *, "Error for sin(x):"
    print *, analytical_solution_sin - arrOut_sin

    print *, "Analytical solution for x^4:"
    print *, analytical_solution_poly
    print *, "numerical solution for x^4:"
    print *, arrOut_poly
    print *, "Error for x^4:"
    print *, analytical_solution_poly - arrOut_poly


    ! Save results to a file
    open(1, file='task2_res.txt', status='replace')
    write(1, *) 'Numerical solution for sin(x):'
    write(1, *) arrOut_sin
    write(1, *) "Analytical solution for sin(x):"
    write(1, *) analytical_solution_sin
    write(1, *) "Error for sin(x):"
    write(1, *) arrOut_sin - analytical_solution_sin

    write(1, *) 'Numerical solution for x^4:'
    write(1, *) arrOut_poly
    write(1, *) "Analytical solution for x^4:"
    write(1, *) analytical_solution_poly
    write(1, *) "Error for x^4:"
    write(1, *) arrOut_poly - analytical_solution_poly
    

    









end program task2









































