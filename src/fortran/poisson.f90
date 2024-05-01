program poisson
    use poisson_solver_utilities
    use csv_writer
    use namelist_utilities
    implicit none
    real, parameter :: pi = 3.14159265358979323846
    integer :: nx=6, ny=6, i, j, iter_count
    double precision :: h, alpha=1, res_rms, f_norm
    double precision, allocatable :: u(:,:), f(:,:)
    character(len=15) :: init_State='spike', output_u='poisson_u', output_f='poisson_f'
    logical :: multigrid=.FALSE.
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
        do j = 1, ny
            f(1, j) = 0.0
            f(nx, j) = 0.0
        end do
        do i = 1, nx
            u(i, 1) = 0.0
            u(i, ny) = 0.0
        end do
    end if
   

    res_rms = sum(abs(f - u)) / (size(f,1) * size(f,2))
    i = 0
    do while (.FALSE.)!(res_rms > 1.0e-5)
        if (multigrid) then
            write(*,*) "Ich bin eine Ente"
            if (mod(i, 50) == 0) then
                if (i < 1000) then
                    write(*,*) "Ich bin eine Ente"
                    call write_u_and_f_to_csv(output_u, output_f, u, f, i)
                end if
            end if
            res_rms = Vcycle_2DPoisson(u,f,h,alpha)
            write(*,*) res_rms  
            i = i + 1
        else
            if (mod(i, 50) == 0) then
                if (i < 1000) then
                    call write_u_and_f_to_csv(output_u, output_f, u, f, i)
                end if
            end if
            res_rms = iteration_2DPoisson(u,f,h,alpha)
            write(*,*) res_rms
            i = i + 1
        end if
    end do

    ! Calculate the norm of f
    f_norm = SQRT(SUM(f**2))

    ! Iterate until residue is less than 1e-5 of f_norm
    iter_count = 0
    res_rms = 1.0D0
    DO WHILE (res_rms> 1.0e-5)
        res_rms = iteration_2DPoisson(u,f,h,alpha)
        iter_count = iter_count + 1
    END DO

    call print_matrix(u)
    call print_matrix(f)
    PRINT *, 'Iterations:', iter_count
    PRINT *, 'Residue:', res_rms

    write(*,*) "Ich bin eine Ente"  ! @Claas das bleibt jetzt drin


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

