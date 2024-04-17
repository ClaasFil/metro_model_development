program HeatDiffusion2d
    use FiniteDifference2d
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none
    integer :: nx, ny, i, j, k, nsteps
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
    outputfilename='heatdif2d_out_random_a03'

    ! Read in Parameters
    ! open(1, file="namelist.nml", status="old")
    ! read(1, *) nx, ny, kappa, total_time, a, outputfilename
    ! close(1)
    ! write(*, *) nx, ny, kappa, total_time, a, outputfilename

    call read_namelist('namelist.nml', nx, ny, kappa, a, outputfilename)

    ! Calculate grid spacing and time step
    h = L / real(nx - 1)
    dt = a * h**2 / kappa
    nsteps = int(total_time / dt)

    ! Allocate memory for temperature arrays
    allocate(T_init(nx, ny), T(nx, ny), T_new(nx, ny), second_derivative(nx, ny))

    ! Initialize temperature array
    T_init = 0.0
    T_init(nx/2, ny/2) = 10.0  ! Spike in the center
    call RANDOM_NUMBER(T_init)
    T = T_init

    ! Write out temperature
    open(unit=10, file=outputfilename, status='unknown')
    write(10, *) T

    ! Integration loop
    do k = 1, int(nsteps)
        ! Compute second derivative
        call finitedifferences(T, h, second_derivative)
        
        ! Update temperature using forward Euler method
        do i = 2, nx - 1
            do j = 2, ny - 1
                T_new(i, j) = T(i, j) + dt * kappa * second_derivative(i, j)
            end do
        end do
        
        ! Apply boundary conditions
        T_new(1, 1) = 0.0
        T_new(nx, ny) = 0.0

        ! Update temperature array for next time step
        T = T_new
        if (mod(k, 10) == 0) then
            ! write(10, *) T
            ! Save temperature data to CSV file
            call write_temperature_to_csv(outputfilename, T, k)
        end if
    end do

    close(10)

    ! Deallocate memory
    deallocate(T_init, T, T_new, second_derivative)

    contains
    subroutine read_namelist(file_path, nx, ny, kappa, a, outfile)
        !! Reads Namelist from given file.
        character(len=*),  intent(in)    :: file_path
        integer,           intent(inout) :: nx, ny
        real,              intent(inout) :: kappa, a
        integer                          :: fu, rc
        character(len=*), intent(inout)  :: outfile

        ! Namelist definition.
        namelist /INPUTS/ nx, ny, kappa, a, outfile

        ! Check whether file exists.
        inquire (file=file_path, iostat=rc)

        if (rc /= 0) then
            write (stderr, '("Error: input file ", a, " does not exist")') file_path
            return
        end if

        ! Open and read Namelist file.
        open (action='read', file=file_path, iostat=rc, newunit=fu)
        read (nml=INPUTS, iostat=rc, unit=fu)
        if (rc /= 0) write (stderr, '("Error: invalid Namelist format")')

        close (fu)
    end subroutine read_namelist

    subroutine write_temperature_to_csv(outfile, T, step)
        real, dimension(:,:), intent(in) :: T
        character(len=*), intent(inout)  :: outfile
        integer, intent(in) :: step
        integer :: i, j
        character(len=100) :: filename
        open(unit=20, file=outfile // trim(adjustl(int2str(step))) // '.csv', status='replace', action='write')
        do i = 1, nx
            do j = 1, ny
                write(20, *) T(i, j)
            end do
        end do
        close(20)
    end subroutine write_temperature_to_csv

    function int2str(i) result(str)
        implicit none
        integer, intent(in) :: i
        character(len=10) :: str

        ! Convert integer to string
        write(str, '(I10)') i
    end function int2str

end program HeatDiffusion2d
