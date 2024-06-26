program nbody
    use nbody_utilities
    use nbody_methods
    implicit none
    integer :: n = 3 ! number of bodies
    integer :: d = 2 ! dimensions
    real :: dt=0.01, t_max=10, t=0 ! time step, maximum time, current time
    integer :: k = 100000 ! number of iterations
    real, parameter :: G = 1.0  ! gravitational constant
    real, allocatable :: x(:, :, :), v(:, :, :), a(:, :, :)  ! positions, velocities and accelerations
    real, allocatable :: m(:) ! masses
    real, allocatable :: xinits(:, :), vinits(:, :) ! initial conditions
    real, allocatable :: x1(:,:), x2(:,:), x3(:,:) ! positions of bodies
    character(len=6) :: inits = "inits1"
    character(len=20) :: method
    integer :: l

    call read_namelist("namelist_nbody.nml", dt, t_max, method, inits)

    ! Read initial conditions from file
    call read_initial_conditions(inits // ".txt", d, n, m, xinits, vinits)

    k = INT(FLOOR(t_max / dt))
    t = 0.0
    
    allocate(x(k, n, d), v(k, n, d), a(k, n, d)) 

    x = 0.0
    v = 0.0
    x(1, :, :) = xinits
    v(1, :, :) = vinits
    
    print *, "Initial positions and velocities"
    call print_matrix(x(1, :, :))
    call print_matrix(v(1, :, :))
    
    print *, "Number of iterations: ", k
    print *, "Method: ", method

    if (method == "FE") then
        call forward_euler(x, v, m, G, dt)
    else if (method == "BE") then
        call backward_euler(x, v, m, G, dt)
    else if (method == "RK4") then
        call runge_kutta4(x, v, m, G, dt)
    end if

    do l = 1, n
        call write_to_csv("output/" // trim(adjustl(method)) // inits // "_body" // trim(adjustl(str(l))) // ".csv", x(:, l, :))
    end do

end program nbody
