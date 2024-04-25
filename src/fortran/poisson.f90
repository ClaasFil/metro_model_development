program HeatDiffusion2d
    implicit none
    real, parameter :: pi = 3.14159265358979323846
    integer :: nx=4, ny=4, i, j, nsteps, k, io_error
    real :: total_time=1., h, dt, kappa=1., a_adv=0.1, a_diff=0.1, B=1.
    real, allocatable :: T_init(:,:), T(:,:), T_new(:,:), d2T(:,:), u(:,:), v(:,:), psi(:,:)
    real, allocatable :: advection(:,:)
    character(len=100) :: init_State='rand', outputfilename, outputfilename_adv
    logical :: multigrid=.TRUE.
    NAMELIST /inputs/ nx, ny, kappa, total_time, a_adv, a_diff, B, init_State

    !read input parameters
    !call read_namelist_ex6('data/namelist/ex_6.nml', nx, ny, kappa, total_time, a_adv, a_diff, B, init_State)
    
    
    ! Init Temperature arrays
    allocate(T_init(nx, ny), T(nx, ny), T_new(nx, ny), d2T(nx, ny), advection(nx, ny))

    write(*,*) "Ich bin ein Zebra"
    
end program HeatDiffusion2d

