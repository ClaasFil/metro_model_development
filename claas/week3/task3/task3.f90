!gfortran -o task3 derivative.f90 task3.f90 && ./task3


program statistics
    use derivative
    implicit none
    integer :: i, L, nx, total_time   !ngth of the model domain (L), the number of grid points (nx), and the total time of the simulation (total_time)
    real :: h, dt                    !grid spacing (h) and the time step (dt)
    real, allocatable :: T_old(:), T_new(:), dT2(:)    !temperature field (T) on the grid


    !stet initial values
    L = 20  ! ength of the model domain
    nx = 20 ! the number of grid points
    total_time = 20 ! the total time of the simulation


    ! calculate the grid spacing
    h = L / (nx - 1)
    dt = 0.4 * h**2

    ! Initialize the temperature on the grid, with a delta function
    allocate(T_old(nx), T_new(nx), dT2(nx))
    dT2 = 0.0
    T_new = 0.0
    T_old = 0.0
    T_old(nx/2) = 1.0

    ! set  the boundary condition T = 0 (for i = 1 and i = nx) independet of initialization
    ! Since the deriverty is calculated only for the inner points, 
    ! the boundary conditions are not updated 
    T_old(1) = 0.0
    T_old(nx) = 0.0


    ! Save results to a file
    open(1, file='task3_res.txt', status='replace')
    write(1, *) 'Temperature field:'
    write(1, *) T_old

    !time loop
    do i = 1, total_time
        
        !calculate the new temperature field
        call deriv2(T_old, dT2, h)
        T_new = T_old + dt * dT2
        write(1, *) T_new
        !print *, T_new
        T_old = T_new

        ! check anf break if every element of T is close to 0
        if (all(abs(T_new) < 0.01)) then
            exit
        end if
    end do


    !deallocate the temperature fields
    deallocate(T_old, T_new, dT2)




end program statistics






