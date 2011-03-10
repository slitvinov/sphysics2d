      subroutine LU_Decomposition(A_input_matrix,b_input_vector,length,
     &                            A_LU_equivalent_matrix,i_MLS_step,
     &                            i_Singular,x_solution,i_in)

      !Solves Ax = b using LU Decomposition
      !A = LU  ==> LUx = b
      !Ux = y  ==> Ly  = b
      !Solve Ly = b using forward subsitution
      !Solve Ux = y using backward substitution
      
      parameter(LU_size_max = 3)
      real A_input_matrix(LU_size_max,LU_size_max)
      real A_temp(LU_size_max,LU_size_max)
      real A_LU_equivalent_matrix(LU_size_max,LU_size_max)
      real L_matrix(LU_size_max,LU_size_max)
      real U_matrix(LU_size_max,LU_size_max)
      real x_solution(LU_size_max), y_solution(LU_size_max)
      real b_input_vector(LU_size_max), sum_Uijxj, sum_Lijyj
      real al_min,al_max
      TOL_REAL = 1.0e-8
      TOL_DBLE = 1.0e-16

      if(i_MLS_step.eq.1)then  !Construct LU Matrix

       !- Check singularity -
       amax_Aij = 0.0
       row_sum_min = 10.0e8
       i_Singular = 0
       do iii = 1,length
         row_sum = 0.0
         do jjj=1,length
           !Transfer input matrix into a temporary matrix
           A_temp(iii,jjj) = A_input_matrix(iii,jjj)
           !Check matrix is not singular
           A_temp_ij_abs = abs(A_temp(iii,jjj))
           row_sum = row_sum + A_temp_ij_abs
           if(A_temp_ij_abs.gt.amax_Aij)then
             amax_Aij = A_temp_ij_abs
           endif
        enddo     
        if(row_sum.lt.row_sum_min)then
          row_sum_min = row_sum
        endif

       enddo
       if(amax_Aij   .lt.TOL_REAL) i_Singular = 1
       if(row_sum_min.lt.TOL_DBLE) i_Singular = 1
       
       if(i_Singular.lt.1)then
         do iii = 1,length
           L_matrix(iii,1)=A_input_matrix(iii,1)
         enddo
         do jjj = 1,length
           U_matrix(1,jjj)=A_input_matrix(1,jjj) / L_matrix(1,1)
         enddo
         !- Forward elimination using Gaussian Elimination for matrices L and U -
         do jjj = 2,length
           do iii = jjj,length
             sum_temp = 0.0
             do kkk = 1,jjj-1
               sum_temp = sum_temp + L_matrix(iii,kkk)*U_matrix(kkk,jjj)
             enddo
             L_matrix(iii,jjj) = A_input_matrix(iii,jjj) - sum_temp
           enddo
           U_matrix(jjj,jjj) = 1.0
           do iii = jjj+1,length
             sum_temp= 0.0
             do kkk = 1,jjj-1
               sum_temp= sum_temp+ L_matrix(jjj,kkk)*U_matrix(kkk,iii)
             enddo
             U_matrix(jjj,iii) = (A_input_matrix(jjj,iii) - sum_temp)
     &                           /L_matrix(jjj,jjj)
           enddo
         enddo 
         !-- Suggestion by Arno Mayrhofer for singular ill-conditioned L matrices --
         al_min = 1.0e8
         al_max = -1.0e8
         do iii=1,length
            if(L_matrix(iii,iii).gt.al_max) al_max = L_matrix(iii,iii)
            if(L_matrix(iii,iii).lt.al_min) al_min = L_matrix(iii,iii)
         end do
         if(abs(al_min).lt.1E-14*abs(al_max)) i_Singular = 1
         !--  End of Arno Suggestion -- 
 
         
         !- Combine into one matrix -
         do iii = 1,length  !Row
           do jjj = iii+1,length !Column
               A_LU_equivalent_matrix(iii,jjj) = U_matrix(iii,jjj)
           enddo
           do jjj = 1,iii !Column
               A_LU_equivalent_matrix(iii,jjj) = L_matrix(iii,jjj)
           enddo
         enddo
       
       endif
  
      else  !if(i_MLS_step.eq.2)then  !Solver equation

      x_solution = 0.0  !zero initialise solution vector
      y_solution = 0.0  !zero initialise temp     vector
      !- Forward Substitution using Lower Matrix -
      y_solution(1)=b_input_vector(1)/A_LU_equivalent_matrix(1,1) !First equation has only one unknown
      do iii = 2,length
        sum_Lijyj = 0.0
        !- Generating the summation term -
        do jjj = 1,iii-1
          sum_Lijyj = sum_Lijyj
     &              + A_LU_equivalent_matrix(iii,jjj)*y_solution(jjj)
        enddo
        y_solution(iii)=(b_input_vector(iii)- sum_Lijyj) /
     &                      A_LU_equivalent_matrix(iii,iii)
      enddo
  
      !- Back Substitution -
      x_solution(length) = y_solution(length)!Nth equation has only one unknown
      do iii = length-1,1,-1
        sum_Uijxj = 0.0
        !- Generating the summation term -
        do jjj = iii+1,length
          sum_Uijxj = sum_Uijxj
     &              + A_LU_equivalent_matrix(iii,jjj)*x_solution(jjj)
        enddo
        x_solution(iii) = (y_solution(iii) - sum_Uijxj)
      enddo
       
      endif   !End of:        if(i_MLS_step.eq.1)then

      end subroutine
