RECURSIVE FUNCTION Vcycle_2DPoisson(u,f,h,alpha) RESULT (res_rms)
use poisson_solver_utilities  !Claas change
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