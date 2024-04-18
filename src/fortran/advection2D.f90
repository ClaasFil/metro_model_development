program HeatDiffusion2d
    ! Import modules
    use constants_module
    use finitedifference
    use namelist_utilities
    use csv_writer
    use matrix_utilities
    use boundaries_ex5
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx=4, ny=4, i, j, nsteps, k, io_error
    real :: total_time=1., h, dt, kappa=1., a_adv=0.1, a_diff=0.1, B=1.
    real, allocatable :: T_init(:,:), T(:,:), T_new(:,:), d2T(:,:), u(:,:), v(:,:), psi(:,:)
    real, allocatable :: advection(:,:)
    character(len=100) :: init_State='rand', outputfilename

    !read input parameters
    call read_namelist_ex5('data/namelist/ex_5.nml', nx, ny, kappa, total_time, a_adv, a_diff, B, init_State)
    
    
    ! Init Temperature arrays
    allocate(T_init(nx, ny), T(nx, ny), T_new(nx, ny), d2T(nx, ny), advection(nx, ny))

    ! Initialize temperature array
    T_init = 0.0
    T = 0.0
    T_new = 0.0
    d2T = 0.0
    advection = 0.0
    
    T_init(nx/2, ny/2) = 10.0  ! Spike in the center
    call RANDOM_NUMBER(T_init)
    T = T_init

   call boundaries_T(T)

   !call print_matrix(T)

    h = 1.0 / real(ny-1)

    ! init velosophi
    allocate(u(nx, ny), v(nx, ny), psi(nx, ny))

    ! Compute the psi matrix
    do i = 1, nx
        do j = 1, ny
            psi(i, j) = B * sin(pi * real(i) / nx) * sin(pi * real(j) / ny)
        end do
    end do

    outputfilename = 'data/ex_5/psi_out.csv'
    open(unit=10, file=trim(outputfilename), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    call write_to_csv(outputfilename, psi)
    close(10)

    

    ! compute u and v
    do i = 2, nx -1
        do j = 2, ny-1
            u(i, j) = (psi(i, j+1) - psi(i, j-1))/2*h
            v(i, j) = (psi(i+1, j) - psi(i-1, j))/2*h
        end do
    end do

    outputfilename = 'data/ex_5/u_out.csv'
    open(unit=10, file=trim(outputfilename), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    call write_to_csv(outputfilename, u)
    close(10)

    outputfilename = 'data/ex_5/v_out.csv'
    open(unit=10, file=trim(outputfilename), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif
    call write_to_csv(outputfilename, v)
    close(10)
    



    !call print_matrix(u)
    !call print_matrix(v)

    ! calc dt dt=MIN(a_diff*h**2/kappa, a_adv*h/vmax)
    dt = MIN(a_diff*h**2/kappa, a_adv*h/MAXVAL(ABS(v)))
    nsteps = int(total_time / dt)

    !call print_matrix(T_init)

    outputfilename = 'data/ex_5/T_out.csv'
    open(unit=10, file=trim(outputfilename), status='replace', action='write', iostat=io_error)
    if (io_error /= 0) then
        print *, 'Error opening file:', io_error
        stop
    endif

    close(10)
    call write_to_csv(outputfilename, T)

    ! Integration loop
    do k = 1, int(nsteps)
        ! Compute second derivative
        call FiniteDifference2D(T, h, d2T)
        call advection2D(T, h, u, v,  advection)
        !call print_matrix(T)
        
        !call print_matrix(d2T)
        !call print_matrix(advection)
        ! apply preconditions to advection 
        !call boundaries_dT(advection)

        !call print_matrix(advection)
        
        ! Update temperature using forward Euler method
        T_new = T + dt * (kappa * d2T - advection)
        call boundaries_T(T_new)
        

        ! Update temperature array for next time step
        T = T_new
        call write_to_csv(outputfilename, T)

        !call print_matrix(T)


        !print *, '------------------------------------------'
        ! Check if 5 steps have been completed
        if (k > 10000) then
            exit
        end if
    end do





    
end program HeatDiffusion2d











