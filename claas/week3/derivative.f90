module derivative
    implicit none
contains
    !Subrutine calculating the secound derivative of an array
    !Input: arrOriginal - array of values
    !       spacing - spacing between the values
    !       derivative2 - random array with the same size as arrOriginal
    !                     where the 2nd derivative will be stored

    subroutine deriv2(arrOriginal, derivative2, spacing)
        implicit none
        real, intent(in) :: spacing, arrOriginal(:)
        real, allocatable, intent(out) :: derivative2(:)
        integer ::n, i
        real, allocatable :: mtx(:,:)

        n = size(arrOriginal)


    ! --------- creat mtx for the 2nd derivative ---------
    ! calculate the 2nd derivative via matrix
    ! creating the matrix if quit expencive but one time cost
    ! calculating the 2nd derivative is cheap and usable many times
    
        allocate(mtx(n, n))

        ! Initialize the matrix with 0s
        mtx = 0.0

        ! Set -2 on the diagonal and 1 above and below the diagonal
        do i = 2, n -1 
            mtx(i,i) = -2.0    ! Set -2 on the diagonal
            if (i < n) then
                mtx(i, i +1) = 1.0 ! Set 1 above the diagonal
            end if
            if (i > 1) then
                mtx(i, i - 1) = 1.0 ! Set 1 above the diagonal
            end if
        end do

        !allocate(derivative2(n))

        derivative2 = matmul(mtx, arrOriginal)/spacing**2
        

    end subroutine deriv2
end module derivative
















