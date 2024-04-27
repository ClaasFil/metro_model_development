program poisson
    use poisson_solver_utilities
    use csv_writer
    use namelist_utilities
    implicit none
    real, parameter :: pi = 3.14159265358979323846
    integer :: nx=16, ny=16, io_error
    real :: h, alpha=0.1, res_rms
    real, allocatable :: u(:,:), f(:,:)
    character(len=100) :: init_State='rand', output_u='poisson_u.csv', output_f='poisson_f.csv'
    logical :: multigrid=.TRUE.
    NAMELIST /inputs/ nx, ny, init_State, multigrid, alpha, output_u, output_f

    !read input parameters
    !call read_namelist_ex6('data/namelist/ex_6.nml', nx, ny, init_State, multigrid, alpha, output_u, output_f)
    
    ! Init variables
    h = 1./real(ny-1)
    
    ! Init arrays
    allocate(u(nx, ny), f(nx, ny))

    u = 0.0

    ! Initialize u
    f = 0.0
    f(nx/2, ny/2) = 10.0  ! Spike in the center
    if (init_State == 'rand') then
        call RANDOM_NUMBER(f)
    end if

    res_rms = sum(abs(f - u)) / (size(f,1) * size(f,2))
    do while (res_rms > 1.0e-5)
        if (multigrid) then
            res_rms = Vcycle_2DPoisson(u,f,h,alpha)
        else
            res_rms = iteration_2DPoisson(u,f,h,alpha)
        end if
    end do

    call write_to_csv(output_u, u)
    call write_to_csv(output_f, f)

    write(*,*) "Ich bin ein Zebra"
    
end program poisson

