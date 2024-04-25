program statistics
    use derivative2D
    use matrix_utilities
    use temperature_io
    use norm_calculations

    implicit none
    INTEGER :: nx=10, ny=10, nsteps, temp=0, iter_count
    REAL :: kappa=0.1, total_time=1.0, a=1.0, h, dt, norm
    real, allocatable :: T_old(:, :), T_new(:, :), dT2(:, :), T_rand(:,:)    !temperature field (T) on the grid
    character(len=100) :: filename  
    

    
    





    ! calculate the grid spacing
    h = 1. / (ny - 1)
    dt = a * h**2 / (kappa)
    nsteps = int(total_time / dt)



    allocate(T_old(nx,ny), T_new(nx,ny), dT2(nx,ny), T_rand(nx,ny))
    dT2 = 0.0
    T_new = 0.0
    T_old = 0.0
    T_old(int(nx/2), int(ny/2)) = 10.0   ! maybe inprofe finding the middle of the grid
    call RANDOM_NUMBER(T_rand)


    
    !call print_matrix(T_old)
   

    ! set  the boundary condition T = 0 (for i = 1 and i = nx) independet of initialization
    ! Since the deriverty is calculated only for the inner points, 
    ! the boundary conditions are not updated 
    T_old(1,:) = 0.0
    T_old(nx,:) = 0.0
    T_old(:,1) = 0.0
    T_old(:,ny) = 0.0
    T_rand(1,:) = 0.0
    T_rand(nx,:) = 0.0
    T_rand(:,1) = 0.0
    T_rand(:,ny) = 0.0

    ! Construct the filename
    !write(filename, '(A, F8.4, A)') 'a_data/', a, '.csv'

    ! Save results to a file
    
    open(1, file='temp.csv', status='replace')
    
    write(1, *) 'norm to a certain a:'
    

    norm = 0.0
    call calc_linfty_norm(T_old, norm)
    write(1, *) 1, ',', norm
    



    close(1)
    deallocate(T_old, T_new, dT2, T_rand)

end program statistics


