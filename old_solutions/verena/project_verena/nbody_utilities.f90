
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


    subroutine write_to_csv(outfile, T)
        implicit none
        real, dimension(:,:), intent(in) :: T
        character(len=*), intent(in)  :: outfile
        integer :: io, i, j

        ! Open the file for appending; create a new one if it doesn't exist
        open(unit=20, file=trim(outfile), status='unknown', action='write', iostat=io, position='append')
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
        open(unit=20, file=trim(outfile), status='unknown', action='write', iostat=io, position='append')
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


    subroutine read_namelist(file_path, nx, ny, kappa, a, outfile)
        use, intrinsic :: iso_fortran_env, only: stderr => error_unit
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

end module nbody_utilities
