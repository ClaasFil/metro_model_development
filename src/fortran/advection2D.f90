program HeatDiffusion2d
    ! Import modules
    use constants_module
    use finitedifference
    use namelist_utilities
    use csv_writer
    use matrix_utilities
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx, ny, i, j, nsteps
    real :: total_time, h, dt, kappa, a_adv, a_diff, B
    real, allocatable :: T_init(:,:), T(:,:), T_new(:,:), second_derivative(:,:), u(:,:), v(:,:), psi(:,:)
    character(len=100) :: outputfilename

    !call read_namelist('data/namelist/namelist_ex5.nml', nx, ny, kappa, total_time, a_adv, a_diff, B, outputfilename)
    ! Initialise Parameters
    nx = 5
    ny = 5
    kappa = 1.0
    total_time = 0.1
    a_adv = 0.03
    a_diff = 0.03
    B = 1
    outputfilename='data/ex5.csv'


    
    ! Init Temperature arrays
    allocate(T_init(nx, ny), T(nx, ny), T_new(nx, ny), second_derivative(nx, ny))

    ! Initialize temperature array
    T_init = 0.0
    T = 0.0
    T_new = 0.0
    second_derivative = 0.0
    
    T_init(nx/2, ny/2) = 10.0  ! Spike in the center
    call RANDOM_NUMBER(T_init)
    !T = T_init

    h = 1.0 / real(ny-1)

    ! init velosophi
    allocate(u(nx, ny), v(nx, ny), psi(nx, ny))

    ! Compute the psi matrix
    do i = 1, nx
        do j = 1, ny
            psi(i, j) = B * sin(pi * real(i) / nx) * sin(pi * real(j) / ny)
        end do
    end do

    call print_matrix(psi)


    ! compute u and v
    do i = 2, nx -1
        do j = 2, ny-1
            u(i, j) = (psi(i, j+1) - psi(i, j-1))/2*h
            v(i, j) = (psi(i+1, j) - psi(i-1, j))/2*h
        end do
    end do

    call print_matrix(u)
    call print_matrix(v)

    ! calc dt dt=MIN(a_diff*h**2/kappa, a_adv*h/vmax)
    dt = MIN(a_diff*h**2/kappa, a_adv*h/MAXVAL(ABS(v)))
    nsteps = int(total_time / dt)





    
end program HeatDiffusion2d











