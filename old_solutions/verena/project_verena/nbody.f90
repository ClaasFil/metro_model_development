program nbody
    use nbody_utilities
    implicit none
    integer, parameter :: n = 3 ! number of bodies
    integer, parameter :: d = 2 ! dimensions
    real :: dt, t_max, t
    integer, parameter :: k = 1000 ! number of iterations
    real, parameter :: G = 1  ! gravitational constant
    real, allocatable :: x(:, :, :), v(:, :, :), a(:, :, :)  ! positions, velocities and accelerations
    real, allocatable :: m(:) ! masses
    real, allocatable :: xinits(:, :), vinits(:, :), inits(:, :, :) ! initial conditions
    real, allocatable :: x1(:), x2(:), x3(:) ! positions of bodies
    integer :: i, j, l

    allocate(x(k, n, d), v(k, n, d), a(k, n, d), m(n), xinits(n, d), vinits(n, d), inits(n, d, 2))

    ! Set up initial conditions
    xinits(1, 1) = -0.602885898116520
    xinits(1, 2) = 1.059162128863347 - 1
    xinits(2, 1) = 0.252709795391000
    xinits(2, 2) = 1.058254872224370 - 1
    xinits(3, 1) = -0.355389016941814
    xinits(3, 2) = 1.038323764315145 - 1
    vinits(1, 1) = 0.122913546623784
    vinits(1, 2) = 0.747443868604908
    vinits(2, 1) = -0.019325586404545
    vinits(2, 2) = 1.369241993562101
    vinits(3, 1) = -0.103587960218793
    vinits(3, 2) = -2.116685862168820
    x(1, :, :) = xinits
    v(1, :, :) = vinits
    call print_matrix(x(1, :, :))
    call print_matrix(v(1, :, :))

    m = 1.0
    dt = 0.01
    t_max = dt*k
    t = 0.0
    l = 0

    call forward_euler(x, v, a, m, G, dt)

    
    contains

    subroutine calculate_acceleration(x, m, G, a)
        real, intent(in) :: x(:, :), G
        real, intent(in) :: m(:) ! masses
        real, intent(out) :: a(:, :) ! accelerations
        integer :: n, d, i, j, k

        n = size(x, 1)
        d = size(x, 2)
        a = 0.0
        do i = 1, n
            do j = 1, n
                if (i /= j) then
                    do k = 1, d
                        a(i, k) = a(i, k) - G * m(j) * (x(i, k) - x(j, k)) / (norm2(x(i,:) - x(j,:)))**(3)
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
        
        k = size(x, 1)
        n = size(x, 2)
        d = size(x, 3)

        !call write_vector_to_csv("fe_body1.csv", x(1, 1, :))
        !call write_vector_to_csv("fe_body2.csv", x(1, 2, :))
        !call write_vector_to_csv("fe_body3.csv", x(1, 3, :))

        do i=1, k ! iterations
            do j=1, n ! bodies
                call calculate_acceleration(x(i, :, :), m, G, a(i, :, :))
                !print *, "Body ", j, " at iteration ", i, " with position, velocity and acceleration:"
                !print *, "a ", a(i, j, :), "\n"
                !print *, "v ", v(i, j, :), "\n"
                !print *, "x ", x(i, j, :), "\n"
                x(i+1, j, :) = x(i, j, :) + v(i, j, :) * dt + 0.5 * a(i, j, :) * dt**2
                v(i+1, j, :) = v(i, j, :) + a(i, j, :) * dt
            end do
            x1 = x(i+1, 1, :)
            x2 = x(i+1, 2, :)
            x3 = x(i+1, 3, :)
            call write_vector_to_csv("fe_body1.csv", x1)
            call write_vector_to_csv("fe_body2.csv", x2)
            call write_vector_to_csv("fe_body3.csv", x3)
        end do
        
    end subroutine forward_euler

end program nbody
