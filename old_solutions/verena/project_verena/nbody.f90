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
    real, allocatable :: xinits(:, :), vinits(:, :), inits(:, :, :) ! initial conditions
    real, allocatable :: x1(:,:), x2(:,:), x3(:,:) ! positions of bodies
    real :: xinits11=-1.0, xinits12=0.0, xinits21=1.0, xinits22=0.0, xinits31=0.0, xinits32=0.0
    real :: vinits11=0.347111, vinits12=0.532728, vinits21=0.347111, vinits22=0.532728, vinits31=-0.694222, vinits32=-1.065456

    call read_namelist("namelist_nbody1.nml", d, n, dt, t_max, xinits11, xinits12, xinits21, xinits22, xinits31, xinits32, &
    vinits11, vinits12, vinits21, vinits22, vinits31, vinits32)

    k = t_max / dt
    t = 0.0
    
    allocate(x(k, n, d), v(k, n, d), a(k, n, d), m(n), xinits(n, d), vinits(n, d), inits(n, d, 2))

    ! Set up initial conditions

    xinits(1, 1) = xinits11
    xinits(1, 2) = xinits12
    xinits(2, 1) = xinits21
    xinits(2, 2) = xinits22
    xinits(3, 1) = xinits31
    xinits(3, 2) = xinits32
    vinits(1, 1) = vinits11
    vinits(1, 2) = vinits12
    vinits(2, 1) = vinits21
    vinits(2, 2) = vinits22
    vinits(3, 1) = vinits31
    vinits(3, 2) = vinits32

    x = 0.0
    v = 0.0
    x(1, :, :) = xinits
    v(1, :, :) = vinits
    print *, "Initial positions and velocities"
    call print_matrix(x(1, :, :))
    call print_matrix(v(1, :, :))

    m = 1.0
    
    print *, "Number of iterations: ", k

    call forward_euler(x, v, a, m, G, dt)

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

    subroutine forward_euler(x, v, a, m, G, dt)
        real, intent(inout) :: x(:, :, :), v(:, :, :), a(:, :, :) ! positions, velocities, accelerations and masses
        real, intent(in) :: G, dt, m(:) ! gravitational constant, time step and masses
        integer :: k ! number of iterations
        integer :: n ! number of bodies
        integer :: d ! dimensions
        integer :: i, j
        
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)

        do i=1, k ! iterations
            call calculate_acceleration(x(i, :, :), m, G, a(i, :, :))
            do j=1, n ! bodies
                x(i+1, j, :) = x(i, j, :) + v(i, j, :) * dt + 0.5 * a(i, j, :) * dt**2
                v(i+1, j, :) = v(i, j, :) + a(i, j, :) * dt
            end do
        end do
        
    end subroutine forward_euler

    subroutine runge_kutta(x, v, a, m, G, dt)
        real, intent(inout) :: x(:, :, :), v(:, :, :), a(:, :, :) ! positions, velocities, accelerations and masses
        real, allocatable :: v1(:, :), v2(:, :) ! intermediate velocities
        real, allocatable :: k1(:, :), k2(:, :) ! intermediate positions
        real, allocatable :: a1(:, :), a2(:, :), a3(:, :) ! intermediate accelerations
        real, intent(in) :: G, dt, m(:) ! gravitational constant, time step and masses
        integer :: k ! number of iterations
        integer :: n ! number of bodies
        integer :: d ! dimensions
        integer :: i, j
        
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)

        allocate(v1(n, d), v2(n, d))
        allocate(k1(n, d), k2(n, d))
        allocate(a1(n, d), a2(n, d), a3(n, d))

        x1 = 0.0
        x2 = 0.0
        a1 = 0.0
        a2 = 0.0
        a3 = 0.0

        do i=1, k ! iterations
            call calculate_acceleration(x(i, :, :), m, G, a1(:, :))
            x1 = x(i, :, :) + dt/2 * a1(:, :)
            call calculate_acceleration(x1(:, :), m, G, a2(:, :))
            x2 = x(i, :, :) + dt * (-a1(:, :) + 2*a2(:, :))
            call calculate_acceleration(x2(:, :), m, G, a3(:, :))
            v1 = v(i, :, :) + dt/2 * a1(:, :)
            v2 = v(i, :, :) + dt * (-a1(:, :) + 2*a2(:, :))
            do j=1, n ! bodies
                v(i+1, j, :) = v(i, j, :) + dt/6 * (a1(j, :) + 4*a2(j, :) + a3(j, :))
                x(i+1, j, :) = x(i, j, :) + dt/6 * (v(i, j, :) + 4*v1(j, :) + 2*v2(j, :))
            end do
        end do
        
    end subroutine runge_kutta

end program nbody
