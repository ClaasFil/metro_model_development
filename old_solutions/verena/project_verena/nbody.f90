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

    call runge_kutta4(x, v, m, G, dt)

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

    subroutine runge_kutta2(x, v, m, G, dt)
        real, intent(inout) :: x(:, :, :), v(:, :, :)
        real, allocatable :: x1(:, :), x2(:, :), x3(:, :), x4(:, :)
        real, allocatable :: v1(:, :), v2(:, :), v3(:, :), v4(:, :)
        real, allocatable :: a1(:, :), a2(:, :), a3(:, :), a4(:, :)
        real, intent(in) :: G, dt, m(:)
        integer :: k, n, d
        integer :: i, j
    
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)
    
        allocate(x1(n, d), x2(n, d), x3(n, d), x4(n, d))
        allocate(v1(n, d), v2(n, d), v3(n, d), v4(n, d))
        allocate(a1(n, d), a2(n, d), a3(n, d), a4(n, d))
    
        do i=1, k-1 ! iterations
            call calculate_acceleration(x(i, :, :), m, G, a1)
    
            v1 = v(i, :, :)
            x1 = x(i, :, :) + dt/2 * v1
    
            call calculate_acceleration(x1, m, G, a2)
            v2 = v(i, :, :) + dt/2 * a1
            x2 = x(i, :, :) + dt/2 * v2
    
            call calculate_acceleration(x2, m, G, a3)
            v3 = v(i, :, :) + dt/2 * a2
            x3 = x(i, :, :) + dt/2 * v3
    
            call calculate_acceleration(x3, m, G, a4)
            v4 = v(i, :, :) + dt * a3
            x4 = x(i, :, :) + dt * v3
    
            do j=1, n
                v(i+1, j, :) = v(i, j, :) + dt/6 * (a1(j, :) + 2*a2(j, :) + 2*a3(j, :) + a4(j, :))
                x(i+1, j, :) = x(i, j, :) + dt/6 * (v1(j, :) + 2*v2(j, :) + 2*v3(j, :) + v4(j, :))
            end do
        end do
    end subroutine runge_kutta2
    
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
            ! Step 1
            call calculate_acceleration(x(i,:,:), m, G, a1)
            kx1 = v(i,:,:)
            kv1 = a1
            
            ! Step 2
            call calculate_acceleration(x(i,:,:) + 0.5*dt*kx1, m, G, a2)
            kx2 = v(i,:,:) + 0.5*dt*kv1
            kv2 = a2
            
            ! Step 3
            call calculate_acceleration(x(i,:,:) + 0.5*dt*kx2, m, G, a3)
            kx3 = v(i,:,:) + 0.5*dt*kv2
            kv3 = a3
            
            ! Step 4
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
