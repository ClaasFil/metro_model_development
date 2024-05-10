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
            double precision :: res_rms, res
            double precision, allocatable :: R(:,:)
            integer :: i, j, nx, ny
        
            nx = SIZE(u, 1)
            ny = SIZE(u, 2)
            
            res = 0.0D0
            res_rms = 0.0D0
            
            ! Perform Gauss-Seidel iteration
            DO j = 2, ny - 1 
                DO i = 2, nx - 1

                    res = f(i, j) - 1./h**2 * (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j))
                    u(i, j) = u(i, j) - alpha * res * h**2 / 4.0
                    res_rms = res_rms + res**2

                END DO
            END DO
            
            ! Calculate root mean square residue
            res_rms = SQRT(res_rms / real((nx-2) * (ny-2), KIND=8))
            
            RETURN
        END FUNCTION iteration_2DPoisson


        subroutine residue_2DPoisson(u,f,h,res_f)
            double precision, INTENT(OUT) :: res_f(:,:)
            double precision, INTENT(IN) :: u(:,:), f(:,:), h
            integer :: i, j, nx, ny

            nx = SIZE(u, 1)
            ny = SIZE(u, 2)
            
            res_f = 0.0D0
            
            ! Calcualte residual
            DO j = 2, ny - 1 
                DO i = 2, nx - 1
                    res_f(i, j) = f(i, j) - 1./h**2 * (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1) - 4.0 * u(i, j))
                END DO
            END DO

        end subroutine residue_2DPoisson


        SUBROUTINE restrict(fine,coarse)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: fine(:,:)
            REAL(8), INTENT(OUT) :: coarse(:,:)
            INTEGER :: nx, ny, nxc, nyc, i, j
          
            nx = SIZE(fine,1); ny=SIZE(fine,2)
            nxc = (nx+1)/2; nyc = (ny+1)/2
          
            coarse = 0.
          
            DO CONCURRENT (i=1:nxc, j=1:nyc)
              coarse(i,j) = fine(2*i-1,2*j-1)
            END DO
        END SUBROUTINE restrict


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
                        ! corr_f(2*i,2*j) = (corr_c(i,j) + corr_c(i+1, j+1))/2
                        corr_f(2*i,2*j) = (corr_c(i,j) + corr_c(i+1, j+1) + corr_c(i+1,j) + corr_c(i, j+1))/4
                    end if
                end do
            end do
        end subroutine prolongate


        SUBROUTINE prolongate2(coarse,fine)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: coarse(:,:)
            REAL(8), INTENT(OUT) :: fine(:,:)
            INTEGER :: nx, ny, nxc, nyc, i, j
          
            nxc = SIZE(coarse,1); nyc=SIZE(coarse,2)
            nx = 2*nxc-1; ny = 2*nyc-1
          
            fine = 0.
          
            DO CONCURRENT (i=1:nx, j=1:ny)
              fine(i,j) = (coarse(FLOOR  ((i+1)/2.), FLOOR  ((j+1)/2.)) &
                         + coarse(FLOOR  ((i+1)/2.), CEILING((j+1)/2.)) &
                         + coarse(CEILING((i+1)/2.), FLOOR  ((j+1)/2.)) &
                         + coarse(CEILING((i+1)/2.), CEILING((j+1)/2.))) / 4.
            END DO
          END SUBROUTINE prolongate2


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