program nbody
    use nbody_utilities
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

    x1 = x(:, 1, :)
    x2 = x(:, 2, :)
    x3 = x(:, 3, :)
    call write_to_csv("fe_body1.csv", x1)
    call write_to_csv("fe_body2.csv", x2)
    call write_to_csv("fe_body3.csv", x3)

    
    contains

    subroutine calculate_acceleration(x, m, G, a)
        real, intent(in) :: x(:, :), G
        real, intent(in) :: m(:) ! masses
        real, intent(out) :: a(:, :) ! accelerations
        integer :: n, d, b1, b2, dim

        n = size(x, 1)
        d = size(x, 2)
        a = 0.0
        do b1 = 1, n
            do b2 = 1, n
                if (b1 /= b2) then
                    do dim = 1, d
                        a(b1, dim) = a(b1, dim) - G * m(b2) * (x(b1, dim) - x(b2, dim)) / (norm2(x(b1,:) - x(b2,:)))**(3)
                    end do
                end if
            end do
        end do
    end subroutine calculate_acceleration

    subroutine forward_euler(x, v, m, G, dt)
        real, intent(inout) :: x(:, :, :), v(:, :, :) ! positions, velocities, accelerations and masses
        real, intent(in) :: G, dt, m(:) ! gravitational constant, time step and masses
        real, allocatable :: a(:, :, :) ! accelerations
        integer :: k ! number of iterations
        integer :: n ! number of bodies
        integer :: d ! dimensions
        integer :: i, j
        
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)

        allocate(a(k, n, d))

        do i=1, k ! iterations
            call calculate_acceleration(x(i, :, :), m, G, a(i, :, :))
            do j=1, n ! bodies
                x(i+1, j, :) = x(i, j, :) + v(i, j, :) * dt + 0.5 * a(i, j, :) * dt**2
                v(i+1, j, :) = v(i, j, :) + a(i, j, :) * dt
            end do
        end do
        
    end subroutine forward_euler

    subroutine backward_euler(x, v, m, G, dt)
        real, intent(inout) :: x(:, :, :), v(:, :, :) ! positions, velocities
        real, intent(in) :: G, dt, m(:) ! gravitational constant, time step, and masses
        real, allocatable :: a(:, :, :), a_next(:, :) ! accelerations
        real, allocatable :: x_next(:, :), v_next(:, :) ! next positions and velocities
        integer :: k ! number of iterations
        integer :: n ! number of bodies
        integer :: d ! dimensions
        integer :: i, j, iter
        real :: norm_diff
        real, parameter :: tol = 1e-6 ! tolerance for fixed-point iteration
        integer, parameter :: max_iter = 100 ! maximum iterations for fixed-point iteration
    
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)
    
        allocate(a(k, n, d), a_next(n, d), x_next(n, d), v_next(n, d))
    
        do i=1, k-1 ! iterations
            ! Initial guess for the next velocities and positions
            call calculate_acceleration(x(i, :, :), m, G, a(i, :, :))
            x_next = x(i, :, :)
            v_next = v(i, :, :)
    
            ! Fixed-point iteration to solve the implicit equations
            do iter = 1, max_iter
                call calculate_acceleration(x_next, m, G, a_next)
                norm_diff = 0.0
                do j=1, n ! bodies
                    ! Update velocities and positions
                    v_next(j, :) = v(i, j, :) + dt * a_next(j, :)
                    x_next(j, :) = x(i, j, :) + dt * v_next(j, :)
    
                    ! Calculate the norm difference for convergence check
                    norm_diff = norm_diff + sum((v_next(j, :) - v(i, j, :))**2) + sum((x_next(j, :) - x(i, j, :))**2)
                end do
    
                ! Check for convergence
                if (sqrt(norm_diff) < tol) exit
            end do
    
            ! Assign converged values to the arrays
            x(i+1, :, :) = x_next
            v(i+1, :, :) = v_next
            a(i+1, :, :) = a_next
        end do

    end subroutine backward_euler
    
    subroutine runge_kutta4(x, v, m, G, dt)
        implicit none
        real, intent(inout) :: x(:, :, :), v(:, :, :) ! positions, velocities
        real, intent(in) :: m(:), G, dt ! masses, gravitational constant, time step
        integer :: k, n, d ! number of iterations, number of bodies, dimensions
        integer :: i, j ! iteration indices
        real, allocatable :: kx1(:,:), kx2(:,:), kx3(:,:), kx4(:,:)
        real, allocatable :: kv1(:,:), kv2(:,:), kv3(:,:), kv4(:,:)
        real, allocatable :: a1(:,:), a2(:,:), a3(:,:), a4(:,:)
        
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)
    
        allocate(kx1(n,d), kx2(n,d), kx3(n,d), kx4(n,d))
        allocate(kv1(n,d), kv2(n,d), kv3(n,d), kv4(n,d))
        allocate(a1(n,d), a2(n,d), a3(n,d), a4(n,d))
    
        do i=1, k-1 ! iterations
            call calculate_acceleration(x(i,:,:), m, G, a1)
            kx1 = v(i,:,:)
            kv1 = a1
            
            call calculate_acceleration(x(i,:,:) + 0.5*dt*kx1, m, G, a2)
            kx2 = v(i,:,:) + 0.5*dt*kv1
            kv2 = a2

            call calculate_acceleration(x(i,:,:) + 0.5*dt*kx2, m, G, a3)
            kx3 = v(i,:,:) + 0.5*dt*kv2
            kv3 = a3

            call calculate_acceleration(x(i,:,:) + dt*kx3, m, G, a4)
            kx4 = v(i,:,:) + dt*kv3
            kv4 = a4
            
            do j=1, n ! update positions and velocities for each body
                x(i+1,j,:) = x(i,j,:) + dt/6.0 * (kx1(j,:) + 2.0*kx2(j,:) + 2.0*kx3(j,:) + kx4(j,:))
                v(i+1,j,:) = v(i,j,:) + dt/6.0 * (kv1(j,:) + 2.0*kv2(j,:) + 2.0*kv3(j,:) + kv4(j,:))
            end do
        end do

        end subroutine runge_kutta4    

end program nbody
