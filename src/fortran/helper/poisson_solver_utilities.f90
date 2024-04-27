module poisson_solver_utilities
    ! use FiniteDifference
    implicit none

    contains

    RECURSIVE FUNCTION iteration_2DPoisson(u,f,h,alpha) RESULT (res_rms)
        IMPLICIT NONE
        REAL, INTENT(INOUT) :: u(:,:)
        REAL, INTENT(IN) :: f(:,:)
        REAL, INTENT(IN) :: h, alpha
        REAL :: res_rms
        REAL, allocatable :: u_derived(:,:), R(:,:)
    
        allocate(u_derived(size(u,1), size(u,2)), R(size(u,1), size(u,2)))
    
        call FiniteDifference2D(u, h, u_derived)
        R = f - u_derived
    
        u = u - alpha * R * h**2 / 4
    
        res_rms = sqrt(sum(R**2) / size(R,1) / size(R,2))
    
        RETURN
    
    END FUNCTION iteration_2DPoisson


        subroutine residue_2DPoisson(u,f,h,res_f)
            REAL, INTENT(INOUT) :: u(:,:), res_f(:,:)
            REAL, INTENT(IN) :: f(:,:), h
            REAL, allocatable :: u_derived(:,:)

            allocate(u_derived(size(u,1), size(u,2)))

            call FiniteDifference2D(u, h, u_derived)
            res_f = f - u_derived

        end subroutine residue_2DPoisson


        subroutine restrict(res_f,res_c)
            REAL, INTENT(IN) :: res_f(:,:)
            REAL, allocatable, INTENT(OUT) :: res_c(:,:)
            INTEGER :: i, j

            allocate(res_c(size(res_f,1)/2, size(res_f,2)/2))

            res_c = 0
            do j = 2, size(res_f,2)-1, 2
                do i = 2, size(res_f,1)-1, 2
                    res_c(i/2,j/2) = res_f(i,j) !+ res_f(i+1,j) + res_f(i,j+1) + res_f(i+1,j+1)
                end do
            end do
        end subroutine restrict


        subroutine prolongate(corr_c,corr_f)
            REAL, INTENT(IN) :: corr_c(:,:)
            REAL, allocatable, INTENT(OUT) :: corr_f(:,:)
            INTEGER :: i, j

            allocate(corr_f(2*size(corr_c,1)-1, 2*size(corr_c,2)-1))
        
            corr_f = 0
            do j = 1, size(corr_c,2)
                do i = 1, size(corr_c,1)
                    corr_f(2*i-1,2*j-1) = corr_c(i,j)
                    if (i < size(corr_c,1)) then
                        corr_f(2*i,2*j-1) = (corr_c(i,j) + corr_c(i+1, j))/2
                    end if
                    if (j < size(corr_c,2)) then
                        corr_f(2*i-1,2*j) = (corr_c(i,j) + corr_c(i, j+1))/2
                    end if
                    if (i < size(corr_c,1) .and. j < size(corr_c,2)) then
                        corr_f(2*i,2*j) = (corr_c(i,j) + corr_c(i+1, j+1))/2
                    end if
                end do
            end do
        end subroutine prolongate


        RECURSIVE FUNCTION Vcycle_2DPoisson(u,f,h,alpha) RESULT (res_rms)
            IMPLICIT NONE
            REAL, INTENT(INOUT) :: u(:,:)
            REAL, INTENT(IN) :: f(:,:), h, alpha
            REAL :: res_rms ! root mean square residue
            INTEGER :: nx, ny, nxc, nyc, i ! local variables
            REAL, ALLOCATABLE :: res_c(:,:), corr_c(:,:), res_f(:,:), corr_f(:,:)
            
            nx = SIZE(u,1); ny = SIZE(u,2) ! must be power of 2 plus 1
            nxc = (nx+1)/2; nyc = (ny+1)/2 ! coarse grid
            
            IF (MIN(nx,ny) > 5) THEN
            
            ALLOCATE(res_f(nx,ny), corr_f(nx,ny))
            ALLOCATE(res_c(nxc,nyc), corr_c(nxc,nyc))
            
            ! take two iterations on the fine grid
            res_rms = iteration_2DPoisson(u,f,h,alpha)
            res_rms = iteration_2DPoisson(u,f,h,alpha)
            
            ! restrict residue to the coarse grid
            CALL residue_2DPoisson(u,f,h,res_f)
            CALL restrict(res_f,res_c)
            
            ! solve for the coarse grid correction
            corr_c = 0.
            res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2,alpha) ! recursion
            
            ! prolongate (interpolate) the correction to the fine grid
            CALL prolongate(corr_c,corr_f)
            
            ! correct the fine-grid solution
            u = u + corr_f
            
            ! take two more smoothing iterations on the fine grid
            res_rms = iteration_2DPoisson(u,f,h,alpha)
            res_rms = iteration_2DPoisson(u,f,h,alpha)
            
            DEALLOCATE(res_f,corr_f,res_c,corr_c)
            
            ELSE ! coarsest level (ny=5): iterate to get "exact" solution
            
            DO i = 1,100
                res_rms = iteration_2DPoisson(u,f,h,alpha)
            END DO
            
            END IF
        
        END FUNCTION Vcycle_2DPoisson


        subroutine FiniteDifference2D(y, h, second_derivative)
            real, intent(in) :: y(:,:)
            real, intent(in) :: h
            real, intent(out) :: second_derivative(:,:)
            integer :: n, m, i, j
            
            n = size(y, 1)
            m = size(y, 2)
            
    
            second_derivative = 0.0
            
    
            do j = 2, m - 1
                do i = 2, n - 1
                    second_derivative(i,j) = (y(i+1,j) + y(i-1,j) + y(i,j+1) - 4.0 * y(i,j) + y(i,j-1)) / h**2
                end do
            end do
    
            ! Apply boundary conditions
            ! Top and bottom rows
            do j = 1, m
                second_derivative(1,j) = (y(2,j) - y(1,j)) / h**2  ! Top row
                second_derivative(n,j) = (y(n-1,j) - y(n,j)) / h**2  ! Bottom row
            end do
    
            ! Left and right columns
            do i = 2, n - 1
                second_derivative(i,1) = (y(i,2) - y(i,1)) / h**2  ! Left column
                second_derivative(i,m) = (y(i,m-1) - y(i,m)) / h**2  ! Right column
            end do
    
        end subroutine FiniteDifference2D
    
end module poisson_solver_utilities