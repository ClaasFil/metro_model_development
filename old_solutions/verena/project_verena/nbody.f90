program nbody
    implicit none
    integer, parameter :: n = 3
    real :: dt, t_max, t
    real, parameter :: G = 6.67430e-11  ! gravitational constant
    real :: x(n), y(n), vx(n), vy(n), m(n)
    integer :: i, j
    real :: ax(n), ay(n)

    ! Set up initial conditions
    x = [0.0, 1.0, 2.0]
    y = [0.0, 0.0, 0.0]
    vx = [0.0, 0.0, 0.0]
    vy = [0.0, 0.5, -0.5]
    m = [1.0, 1.0, 1.0]
    dt = 0.01
    t_max = 10.0

    ! Perform simulation
    do while (t < t_max)
        ! Calculate accelerations
        do i = 1, n
            ax(i) = 0.0
            ay(i) = 0.0
            do j = 1, n
                if (i /= j) then
                    ax(i) = ax(i) - G * m(j) * (x(i) - x(j)) / ((x(i) - x(j))**2 + (y(i) - y(j))**2)**(3./2)
                    ay(i) = ay(i) - G * m(j) * (y(i) - y(j)) / ((x(i) - x(j))**2 + (y(i) - y(j))**2)**(3./2)
                end if
            end do
        end do
        
        ! Update positions and velocities using forward Euler method
        do i = 1, n
            vx(i) = vx(i) + ax(i) * dt
            vy(i) = vy(i) + ay(i) * dt
            x(i) = x(i) + vx(i) * dt
            y(i) = y(i) + vy(i) * dt
        end do

        ! Increment time
        t = t + dt

        ! Output positions
        print *, t, x, y
    end do

end program nbody
