program dry_convection
    ! Import modules
    ! use constants_module
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
    character(len=100) :: init_State='rand', outputfilename, outputfilename_adv

    !read input parameters
    call read_namelist_ex5('data/namelist/ex_5.nml', nx, ny, kappa, total_time, a_adv, a_diff, B, init_State)
    
    
    ! Init Temperature arrays
    allocate(T_init(nx, ny), T(nx, ny), T_new(nx, ny), d2T(nx, ny), advection(nx, ny))


    print *, 'I am a zebra'





    
end program dry_convection











