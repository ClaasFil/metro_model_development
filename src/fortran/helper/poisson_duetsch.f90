MODULE poisson_duetsch

    IMPLICIT NONE
    INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(15,307)
    
    CONTAINS
    
    RECURSIVE FUNCTION Vcycle_2DPoisson_duetsch(u,f,h,alpha) RESULT (res_rms)
      IMPLICIT NONE
      REAL(KIND=DBL), INTENT(INOUT) :: u(:,:)
      REAL(KIND=DBL), INTENT(IN) :: f(:,:), h, alpha
      REAL(KIND=DBL) :: res_rms ! root mean square residual
      INTEGER :: nx, ny, nxc, nyc, i ! local variables
      REAL(KIND=DBL), ALLOCATABLE :: res_c(:,:), corr_c(:,:), res_f(:,:), corr_f(:,:)
    
      nx = SIZE(u,1); ny = SIZE(u,2) ! must be power of 2 plus 1
      nxc = (nx+1)/2; nyc = (ny+1)/2 ! coarse grid
    
      IF (MIN(nx,ny) > 5) THEN
    
        ALLOCATE(res_f(nx,ny), corr_f(nx,ny))
        ALLOCATE(res_c(nxc,nyc), corr_c(nxc,nyc))
    
        ! take two iterations on the fine grid
        res_rms = iteration_2DPoisson(u,f,h,alpha)
        res_rms = iteration_2DPoisson(u,f,h,alpha)
    
        ! restrict residual to the coarse grid
        CALL residue_2DPoisson(u,f,h,res_f)
        CALL restrict(res_f,res_c)
    
        ! solve for the coarse grid correction
        corr_c = 0.
        res_rms = Vcycle_2DPoisson_duetsch(corr_c,res_c,h*2,alpha) ! recursion
    
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
    
    END FUNCTION Vcycle_2DPoisson_duetsch
    
    FUNCTION iteration_2DPoisson(u,f,h,alpha)
      IMPLICIT NONE
      REAL(KIND=DBL), INTENT(INOUT) :: u(:,:)
      REAL(KIND=DBL), INTENT(IN) :: f(:,:), h, alpha
      REAL(KIND=DBL) :: iteration_2DPoisson, res, sum_res2
      INTEGER :: nx, ny, i, j
    
      nx = SIZE(u,1); ny =  SIZE(u,2)
    
      sum_res2 = 0.
      DO j=2,ny-1
        DO i=2,nx-1
          res = f(i,j) - 1./h**2 * (u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) - 4*u(i,j))
          u(i,j) = u(i,j) - alpha*res*h**2/4.
          sum_res2 = sum_res2 + res**2
        END DO
      END DO
    
      iteration_2DPoisson = SQRT(sum_res2/(nx*ny))
    END FUNCTION iteration_2DPoisson
    
    SUBROUTINE residue_2DPoisson(u,f,h,res)
      IMPLICIT NONE
      REAL(KIND=DBL), INTENT(IN) :: u(:,:), f(:,:), h
      REAL(KIND=DBL), INTENT(OUT) :: res(:,:)
      INTEGER :: nx, ny, i, j
    
      nx = SIZE(u,1); ny = SIZE(u,2)
    
      DO CONCURRENT (i=2:nx-1, j=2:ny-1)
        res(i,j) = f(i,j) - 1./h**2 * (u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) - 4*u(i,j))
      END DO
    END SUBROUTINE residue_2DPoisson
    
    SUBROUTINE restrict(fine,coarse)
      IMPLICIT NONE
      REAL(KIND=DBL), INTENT(IN) :: fine(:,:)
      REAL(KIND=DBL), INTENT(OUT) :: coarse(:,:)
      INTEGER :: nx, ny, nxc, nyc, i, j
    
      nx = SIZE(fine,1); ny=SIZE(fine,2)
      nxc = (nx+1)/2; nyc = (ny+1)/2
    
      coarse = 0.
    
      DO CONCURRENT (i=1:nxc, j=1:nyc)
        coarse(i,j) = fine(2*i-1,2*j-1)
      END DO
    END SUBROUTINE restrict
    
    SUBROUTINE prolongate(coarse,fine)
      IMPLICIT NONE
      REAL(KIND=DBL), INTENT(IN) :: coarse(:,:)
      REAL(KIND=DBL), INTENT(OUT) :: fine(:,:)
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
    END SUBROUTINE prolongate
    
    END MODULE poisson_duetsch