
module matrix_utilities
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
            write(*, '(100F8.2)') (matrix(i, j), j = 1, m)
        end do
    end subroutine print_matrix
end module matrix_utilities




module csv_writer
    implicit none
contains
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
                    write(20, '(F6.2)') T(i, j)  ! Last element in the row
                else
                    write(20, '(F6.2, A)', advance='no') T(i, j), ','  ! Elements with comma
                end if
            end do
            write(20, *)  ! Newline for the next row
        end do

        ! Close the file
        close(20)
    end subroutine write_to_csv
end module csv_writer


module namelist_utilities
    implicit none

    contains

    
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

    subroutine read_namelist_ex5(filename, nx, ny, kappa, total_time, a_adv, a_diff, B, init_State)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: nx, ny
        real, intent(out) :: kappa, total_time, a_adv, a_diff, B
        character(len=*), intent(out) :: init_State
    
        namelist /INPUTS/ nx, ny, kappa, total_time, a_adv, a_diff, B, init_State
    
        integer :: io
        open(unit=10, file=filename, status='old', iostat=io)
        if (io /= 0) then
            print *, "Failed to open namelist file:", filename
            return
        end if
    
        read(unit=10, nml=INPUTS, iostat=io)
        if (io /= 0) then
            print *, "Error reading namelist"
        end if
    
        close(unit=10)
    end subroutine read_namelist_ex5
end module namelist_utilities




module constants_module
    implicit none
    real, parameter :: pi = 3.14159265358979323846
end module constants_module




module boundaries_ex5
    !use matrix_utilities
    implicit none
contains
    subroutine boundaries_T(matrix)
        real, dimension(:, :), intent(inout) :: matrix
        integer :: n, m

        n = size(matrix, 1)
        m = size(matrix, 2)
        
        !call print_matrix(matrix)
        ! Apply boundary conditions T = 1 at y = 0 (j=1), T = 0 at y = ymax (j=ny)
        matrix(:, 1) = 1.0
        matrix(:, m) = 0.0
    
        !call print_matrix(matrix)


    end subroutine boundaries_T


    subroutine boundaries_dT(matrix)
        real, dimension(:, :), intent(inout) :: matrix
        integer :: n, m

        n = size(matrix, 1)
        m = size(matrix, 2)
        
        !call print_matrix(matrix)
        ! Apply boundary conditions Tx = 0 at x = 0 (i=1), Tx = 0 at x = xmax (i=nx)
        matrix(1, :) = 0.0
        matrix(n, :) = 0.0
    
        !call print_matrix(matrix)

    end subroutine boundaries_dT

end module boundaries_ex5




