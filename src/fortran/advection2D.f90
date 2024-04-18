program HeatDiffusion2d
    ! Import modules
    use finitedifference
    use namelist_utilities
    use csv_writer
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx, ny, i, j, k, nsteps, io_error
    real :: L, total_time, h, dt, time, a, kappa
    real, allocatable :: T_init(:,:), T(:,:), T_new(:,:), second_derivative(:,:)
    character(len=30) :: outputfilename

    ! Initialise Parameters
    L = 1.0           ! Length of the domain
    nx = 64
    ny = 64
    kappa = 1.0
    total_time = 0.1
    a = 0.3
    outputfilename='data/heatdif2d_out_random_a03'


    call read_namelist('data/namelist/namelist.nml', nx, ny, kappa, a, outputfilename)


    ! Print a header
    print *, 'List of Zoo Animals:'
    print *,  'Elephant ', 'Giraffe ', 'Lion ', 'Tiger ', 'Bear ', 'Zebra ', 'Panda ', 'Kangaroo ', 'Rhinoceros ', 'Hippopotamus ' 


    
end program HeatDiffusion2d











