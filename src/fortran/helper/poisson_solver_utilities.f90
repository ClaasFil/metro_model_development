module poisson_solver_utilities
    use matrix_utilities
    implicit none

    contains

        FUNCTION iteration_2DPoisson(u,f,h,alpha) RESULT (res_rms)
            use matrix_utilities
            IMPLICIT NONE
            double precision, INTENT(INOUT) :: u(:,:)
            double precision, INTENT(IN) :: f(:,:)
            double precision, INTENT(IN) :: h, alpha
            double precision :: res_rms
            double precision, allocatable :: R(:,:)
            integer :: i, j, nx, ny
        
            nx = SIZE(u, 1)
            ny = SIZE(u, 2)

            allocate(R(size(u,1), size(u,2)))
            R = 0.0D0
            
            res_rms = 0.0

            !call print_matrix(R)
            !call print_matrix(u)
            !call print_matrix(f)
            
            ! Perform Gauss-Seidel iteration
            DO j = 2, ny - 1 
                DO i = 2, nx - 1
                    u(1,:) = 0.0D0
                    u(nx,:) = 0.0D0
                    u(:,1) = 0.0D0
                    u(:,ny) = 0.0D0
                    !R = f(i, j) - (1.0 / h**2) * (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j))
                    !u(i, j) = u(i, j) - alpha * R(i, j)
                    !res_rms = res_rms + R(i, j)**2

                    R(i,j) = f(i, j) - 0.25D0 * (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j))
                    u(i, j) = u(i, j) - alpha * R(i,j) * h**2 * 0.25D0
                    res_rms = res_rms + R(i,j)**2
                    !call print_matrix(u)
                    !PRINT *, 'R:', u(i, j-1) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j)
                    !PRINT *, 'u:', u(i,j)
                    !PRINT *, 'f:', f(i,j)
                    !PRINT *, 'alpha:', alpha
                    !PRINT *, 'h:', h
                END DO
            END DO
            
            ! Calculate root mean square residue
            res_rms = SQRT(res_rms / real((nx-2) * (ny-2), KIND=8))
            
            RETURN
        END FUNCTION iteration_2DPoisson


        subroutine residue_2DPoisson(u,f,h,res_f)
            double precision, INTENT(INOUT) :: res_f(:,:)
            double precision, INTENT(IN) :: u(:,:), f(:,:), h
            integer :: i, j, nx, ny
            double precision :: res_rms

            nx = SIZE(u, 1)
            ny = SIZE(u, 2)
            
            res_rms = 0.0D0

            !call print_matrix(R)
            !call print_matrix(u)
            !call print_matrix(f)
            
            ! Perform Gauss-Seidel iteration
            DO j = 2, ny - 1 
                DO i = 2, nx - 1
                    res_f(i, j) = f(i, j) - (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j))
                    res_rms = res_rms + res_f(i,j)**2
                    !call print_matrix(u)
                    !PRINT *, 'R:', u(i, j-1) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j)
                    !PRINT *, 'u:', u(i,j)
                    !PRINT *, 'f:', f(i,j)
                END DO
            END DO
            
            ! Calculate root mean square residue
            res_rms = SQRT(res_rms / real((nx-2) * (ny-2), KIND=8))

        end subroutine residue_2DPoisson


        subroutine restrict(res_f,res_c)
            double precision, INTENT(IN) :: res_f(:,:)
            double precision, allocatable, INTENT(OUT) :: res_c(:,:)
            INTEGER :: i, j

            allocate(res_c(size(res_f,1)/2, size(res_f,2)/2))

            res_c = 0.0
            do j = 2, size(res_f,2)-1, 2
                do i = 2, size(res_f,1)-1, 2
                    res_c(i/2,j/2) = res_f(i,j) 
                end do
            end do
        end subroutine restrict


        subroutine prolongate(corr_c,corr_f)
            double precision, INTENT(IN) :: corr_c(:,:)
            double precision, allocatable, INTENT(OUT) :: corr_f(:,:)
            INTEGER :: i, j

            allocate(corr_f(2*size(corr_c,1)-1, 2*size(corr_c,2)-1))
        
            corr_f = 0.0D0
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
            double precision, INTENT(INOUT) :: u(:,:)
            double precision, INTENT(IN) :: f(:,:), h, alpha
            double precision :: res_rms ! root mean square residue
            INTEGER :: nx, ny, nxc, nyc, i ! local variables
            double precision, ALLOCATABLE :: res_c(:,:), corr_c(:,:), res_f(:,:), corr_f(:,:)
            
            nx = SIZE(u,1); ny = SIZE(u,2) ! must be power of 2 plus 1
            nxc = (nx+1)/2; nyc = (ny+1)/2 ! coarse grid
            
            IF (MIN(nx,ny) > 5) THEN
            
                ALLOCATE(res_f(nx,ny), corr_f(nx,ny))
                ALLOCATE(res_c(nxc,nyc), corr_c(nxc,nyc))
                
                ! take two iterations on the fine grid
                res_rms = iteration_2DPoisson(u,f,h,alpha)
                res_rms = iteration_2DPoisson(u,f,h,alpha)

                !PRINT *, 'u:'
                !call print_matrix(u)
                
                ! restrict residue to the coarse grid
                CALL residue_2DPoisson(u,f,h,res_f)
                !PRINT *, 'res_f:'
                !call print_matrix(res_f)
                CALL restrict(res_f,res_c)
                !PRINT *, 'res_c:'
                !call print_matrix(res_c)
                
                ! solve for the coarse grid correction
                corr_c = 0.
                res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2,alpha) ! recursion
                !PRINT *, 'Corr_c:'
                !call print_matrix(corr_c)
                
                ! prolongate (interpolate) the correction to the fine grid
                CALL prolongate(corr_c,corr_f)
                !PRINT *, 'Prolongated f:'
                !call print_matrix(corr_f)
                
                ! correct the fine-grid solution
                u = u + corr_f
                
                ! take two more smoothing iterations on the fine grid
                res_rms = iteration_2DPoisson(u,f,h,alpha)
                !PRINT *, 'res_rms:', res_rms
                res_rms = iteration_2DPoisson(u,f,h,alpha)
                !PRINT *, 'res_rms:', res_rms
                
                DEALLOCATE(res_f,corr_f,res_c,corr_c)
                
            ELSE ! coarsest level (ny=5): iterate to get "exact" solution
            
                DO i = 1,100
                    res_rms = iteration_2DPoisson(u,f,h,alpha)
                    !PRINT *, 'res_rms:', res_rms
                END DO
            
            END IF
        
        END FUNCTION Vcycle_2DPoisson
    
end module poisson_solver_utilities