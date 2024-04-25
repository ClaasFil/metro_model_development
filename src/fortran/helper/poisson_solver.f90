module poisson_solver
    implicit none
    use FiniteDifferences

    contains

        FUNCTION iteration_2DPoisson(u,f,h,alpha)
            REAL, INTENT(INOUT) :: u(:,:)
            REAL, INTENT(IN) :: f(:,:), h, alpha
            REAL, INTENT(OUT) :: rms
            REAL(:,:) :: u_derived(:,:), R(:,:)

            call FiniteDifference2D(u, h, u_derived)
            R = f - u_derived
            u = u - alpha * R * h^2 / 4
            rms = sqrt(sum(R**2) / size(R,1) / size(R,2))

            RETURN

        END iteration_2DPoisson


        subroutine residue_2DPoisson(u,f,h,res_f)
            REAL, INTENT(INOUT) :: u(:,:), res_f(:,:)
            REAL, INTENT(IN) :: f(:,:), h
            REAL(:,:) :: u_derived(:,:)

            call FiniteDifference2D(u, h, u_derived)
            res_f = f - u_derived

        end residue_2DPoisson


        subroutine restrict(res_f,res_c)
            REAL, INTENT(IN) :: res_f(:,:)
            REAL, INTENT(OUT) :: res_c(:,:)
            INTEGER :: i, j

            res_c = 0
            do j = 2, size(res_f,2)-1, 2
                do i = 2, size(res_f,1)-1, 2
                    res_c(i/2,j/2) = res_f(i,j) !+ res_f(i+1,j) + res_f(i,j+1) + res_f(i+1,j+1)
                end do
            end do
        end restrict


        subroutine prolongate(corr_c,corr_f)
            REAL, INTENT(IN) :: corr_c(:,:)
            REAL, INTENT(OUT) :: corr_f(:,:)
            INTEGER :: i, j

            corr_f = 0
            corr_f(1,1) = corr_c(1,1)
            corr_f(size(corr_f,1),1) = corr_c(size(corr_c,1),1)
            corr_f(1,size(corr_f,2)) = corr_c(1,size(corr_c,2))
            corr_f(size(corr_f,1),size(corr_f,2)) = corr_c(size(corr_c,1),size(corr_c,2))
            do j = 1, size(corr_c,2)
                do i = 1, size(corr_c,1)
                    corr_f(2*i-1,2*j-1) = (corr_c(i,j) + corr_c(i-1, j-1))/2
                    corr_f(2*i,2*j-1) = (corr_c(i,j) + corr_c(i, j-1))/2
                    corr_f(2*i-1,2*j) = (corr_c(i,j) + corr_c(i-1, j))/2
                    corr_f(2*i,2*j) = corr_c(i,j)
                end do
            end do
        end prolongate
    
end module poisson_solver