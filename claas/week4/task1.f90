program statistics
    use derivative2D
    use matrix_utilities
    implicit none
    INTEGER :: nx=10, ny=10, nsteps, temp=0, iter_count
    REAL :: kappa=0.1, total_time=1.0, a=1.0, h, dt, i
    real, allocatable :: T_old(:, :), T_new(:, :), dT2(:, :)    !temperature field (T) on the grid
    ! Rean in the namelist
    NAMELIST /inputs/ nx, ny, kappa, total_time, a


    open(10, file='param.txt')
    read(10, inputs)
    close(10)




    ! calculate the grid spacing
    h = 1. / (ny - 1)
    dt = a * h**2 / (kappa)
    nsteps = int(total_time / dt)



    allocate(T_old(nx,ny), T_new(nx,ny), dT2(nx,ny))
    dT2 = 0.0
    T_new = 0.0
    T_old = 0.0
    T_old(int(nx/2), int(ny/2)) = 10.0   ! maybe inprofe finding the middle of the grid


    
    !call print_matrix(T_old)
   

    
    

    ! set  the boundary condition T = 0 (for i = 1 and i = nx) independet of initialization
    ! Since the deriverty is calculated only for the inner points, 
    ! the boundary conditions are not updated 
    T_old(1,:) = 0.0
    T_old(nx,:) = 0.0
    T_old(:,1) = 0.0
    T_old(:,ny) = 0.0

    ! Save results to a file
    open(1, file='task1_res.txt', status='replace')
    write(1, *) 'nx'
    write(1, *) nx
    write(1, *) 'ny'
    write(1, *) ny
    write(1, *) T_old

    !time loop

    do while (iter_count < nsteps)
        
        !calculate the new temperature field
        call deriv2(T_old, h, dT2)
        !call print_matrix(T_old + (dt* kappa * dT2))
        T_new = T_old + dt* kappa * dT2
        call print_matrix(T_new)
        write(1, *) T_new
        !print *, T_new
        T_old = T_new

        ! check anf break if every element of T is close to 0
        if (all(abs(T_new) < 0.01)) then
            exit
        end if

        if (temp > 10) then
            exit
        end if
        temp = temp +  1

        iter_count = iter_count + 1
    end do

    deallocate(T_old, T_new, dT2)



end program statistics