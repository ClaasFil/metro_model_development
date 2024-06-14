
module nbody_utilities
    implicit none
contains

    subroutine print_matrix(matrix)
        real, dimension(:, :), intent(in) :: matrix
        integer :: i, j, n, m

        n = size(matrix, 1)
        m = size(matrix, 2)

        ! Print the matrix in a formatted manner
        print *, 'Matrix:'
        do i = 1, n
            ! Print each row with aligned columns
            write(*, '(100F12.2)') (matrix(i, j), j = 1, m)
        end do
    end subroutine print_matrix


    subroutine read_initial_conditions(filename, n_dimensions, n_bodies, m, x, v)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: n_dimensions, n_bodies
        real, allocatable, intent(out) :: m(:)
        real, allocatable, intent(out) :: x(:, :), v(:, :)
        integer :: i, j
        character(len=100) :: line
        integer :: ios
        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if
    
        ! Read number of dimensions and bodies
        read(10, *) n_dimensions, n_bodies
    
        ! Allocate arrays
        allocate(m(n_bodies))
        allocate(x(n_bodies, n_dimensions))
        allocate(v(n_bodies, n_dimensions))
    
        ! Read masses, positions, and velocities
        do i = 1, n_bodies
            read(10, *) m(i), (x(i, j), j=1, n_dimensions), (v(i, j), j=1, n_dimensions)
        end do
    
        close(10)
    end subroutine read_initial_conditions


    subroutine write_to_csv(outfile, T)
        implicit none
        real, dimension(:,:), intent(in) :: T
        character(len=*), intent(in)  :: outfile
        integer :: io, i, j

        ! Open the file for appending; create a new one if it doesn't exist
        open(unit=20, file=trim(outfile), status='replace', action='write', iostat=io, position='append')
        if (io /= 0) then
            print *, "Failed to open file:", trim(outfile)
            return
        end if

        ! Write the matrix to the file
        do i = 1, size(T, 1)
            do j = 1, size(T, 2)
                if (j == size(T, 2)) then
                    write(20, '(F12.2)') T(i, j)  ! Last element in the row
                else
                    write(20, '(F12.2, A)', advance='no') T(i, j), ','  ! Elements with comma
                end if
            end do
            write(20, *)  ! Newline for the next row
        end do

        ! Close the file
        close(20)
    end subroutine write_to_csv


    subroutine write_vector_to_csv(outfile, T)
        implicit none
        real, dimension(:), intent(in) :: T
        character(len=*), intent(in)  :: outfile
        integer :: io, i, j

        ! Open the file for appending; create a new one if it doesn't exist
        open(unit=20, file=trim(outfile), status='replace', action='write', iostat=io, position='append')
        if (io /= 0) then
            print *, "Failed to open file:", trim(outfile)
            return
        end if

        ! Write the matrix to the file
        do i = 1, size(T, 1)
            if (i == size(T, 1)) then
                write(20, '(F12.2)') T(i)  ! Last element in the row
            else
                write(20, '(F12.2, A)', advance='no') T(i), ','  ! Elements with comma
            end if
        end do

        ! Close the file
        close(20)
    end subroutine write_vector_to_csv


    subroutine read_namelist(file_path, dt, t_max, method, inits)
        use, intrinsic :: iso_fortran_env, only: stderr => error_unit
        !! Reads Namelist from given file.
        character(len=*),  intent(in)    :: file_path
        real,           intent(inout) :: dt, t_max
        character(len=20), intent(inout) :: method
        character(len=6), intent(inout) :: inits
        integer :: rc, fu
        !character(len=*), intent(inout)  :: outfile

        ! Namelist definition.
        namelist /INPUTS/ dt, t_max, method, inits

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

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end module nbody_utilities
