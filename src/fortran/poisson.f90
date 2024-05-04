program poisson
    use poisson_solver_utilities
    use csv_writer
    use namelist_utilities
    implicit none
    real, parameter :: pi = 3.14159265358979323846
    integer :: nx=17, ny=17, i, j, iter_count
    double precision :: h, alpha=1, res_rms, f_norm
    double precision, allocatable :: u(:,:), f(:,:)
    character(len=9) :: init_State='spike', output_u='poisson_u', output_f='poisson_f'
    logical :: multigrid=.TRUE.
    NAMELIST /inputs/ nx, ny, init_State, multigrid, alpha, output_u, output_f

    !read input parameters
    !call read_namelist_ex6('data/namelist/ex_6.nml', nx, ny, init_State, multigrid, alpha, output_u, output_f)
    
    ! Init variables
    h = 1./real(ny-1)
    
    ! Init arrays
    allocate(u(nx, ny), f(nx, ny))

    u = 0.0D0

    f = 0.0D0
    f(nx/2, ny/2) = 1.0D0  ! Spike in the center

    ! Random initial state with satisfied boundary conditions
    if (init_State == 'rand') then
        call RANDOM_NUMBER(f)
        f(1,:) = 0.0D0
        f(nx,:) = 0.0D0
        f(:,1) = 0.0D0
        f(:,ny) = 0.0D0
    end if

    ! Calculate the norm of f
    f_norm = SQRT(SUM(f**2)/(nx*ny))

    ! Loop to call solver until convergence
    res_rms = f_norm
    i = 0

    do while (res_rms/f_norm > 1.0e-5)
        if (multigrid) then

            ! Save not every iteration step as csv
            if (mod(i, 1) == 0) then
                if (i < 10) then
                    call write_u_and_f_to_csv(output_u, output_f, u, f, i)
                end if
            end if
            !write(*,*) i
            res_rms = Vcycle_2DPoisson(u,f,h,alpha)
            ! write(*,*) res_rms  

        else

            ! Save not every iteration step as csv
            if (mod(i, 10) == 0) then
                if (i < 100) then
                    call write_u_and_f_to_csv(output_u, output_f, u, f, i)
                end if
            end if

            res_rms = iteration_2DPoisson(u,f,h,alpha)
            ! write(*,*) res_rms

        end if
        i = i + 1
    end do
    write(*,*) i


    ! Simpler loop for testing purposes

    iter_count = 0
    res_rms = 1.0D0
    DO WHILE (.FALSE.) !i = 0, 1000
        res_rms = Vcycle_2DPoisson(u,f,h,alpha)
        iter_count = iter_count + 1
        write(*,*) res_rms
    END DO

    ! call print_matrix(u)
    ! call print_matrix(f)


    ! Subroutine for saving to csv (as in utilities module)

    contains

    subroutine write_u_and_f_to_csv(outfile_u, outfile_f, u, f, step)
        double precision, dimension(:,:), intent(in) :: u, f
        character(len=*), intent(inout)  :: outfile_u, outfile_f
        integer, intent(in) :: step
        integer :: i, j
        open(unit=20, file=outfile_u // trim(adjustl(int2str(step))) // '.csv', status='replace', action='write')
        open(unit=30, file=outfile_f // trim(adjustl(int2str(step))) // '.csv', status='replace', action='write')
        do j = 1, size(u, 2)
            do i = 1, size(u, 1)
                write(20, *) u(i, j)
                write(30, *) f(i, j)
            end do
        end do
        close(20)
        close(30)
    end subroutine write_u_and_f_to_csv

    function int2str(i) result(str)
        implicit none
        integer, intent(in) :: i
        character(len=10) :: str

        ! Convert integer to string
        write(str, '(I10)') i
    end function int2str
    
end program poisson

