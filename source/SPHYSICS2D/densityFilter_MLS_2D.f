c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr Alejandro Crespo, Dr Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
c
c    This file is part of SPHYSICS.
c
c    SPHYSICS is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    SPHYSICS is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.

      subroutine densityFilter
       
      include 'common.2D'

      real Amatrix_temp(3,3),Amatrix_temp2(3,3),Ainv_check(3,3)
      real A_input_matrix(3,3), A_LU_equivalent_matrix(3,3)
      real b_input_vector(3)  , x_solution(3)
      real A_MLS_inverse(3,3) , col(3)
       

      nstart_MLS = nbp1   !nstart  
c     -- Zeroing Variables --         
      do i = nstart_MLS,np
c        -- Kernel ReNormalisation ARRAYS --
         pVol(i) = pm(i)/rhop(i)
         sum_wab(i)  = 0.0
         rhop_sum(i) = 0.0
         rhop_sum_MLS(i) = 0.0
         A_MLS(i,1,1) = 0.0
         A_MLS(i,1,2) = 0.0
         A_MLS(i,1,3) = 0.0
         A_MLS(i,2,1) = 0.0
         A_MLS(i,2,2) = 0.0
         A_MLS(i,2,3) = 0.0
         A_MLS(i,3,1) = 0.0
         A_MLS(i,3,2) = 0.0
         A_MLS(i,3,3) = 0.0
         !- test matrix -
         !A_input_matrix(1,1) =  1.0
         !A_input_matrix(1,2) = -1.0
         !A_input_matrix(1,3) =  2.0
         !A_input_matrix(2,1) =  3.0
         !A_input_matrix(2,2) =  0.0
         !A_input_matrix(2,3) =  1.0
         !A_input_matrix(3,1) =  1.0
         !A_input_matrix(3,2) =  0.0
         !A_input_matrix(3,3) =  2.0
         !b_input_vector(1) = 12.0
         !b_input_vector(2) = 11.0
         !b_input_vector(3) = 2.0
      enddo
                 
      i_MLS_part = 1
      call ac_MLS    !Perform the correction matrix construction

      !construct LU-decomposition matrix
      i_singularMatrix_count = 0
      length = 3
      do i=nstart_MLS,np
        !temporary matrix 
        !Self contribution
        A_MLS(i,1,1) = A_MLS(i,1,1) + adh*pm(i)/rhop(i)  !self contribution
        A_MLS_11 = A_MLS(i,1,1)
        A_MLS(i,2,1) = A_MLS(i,1,2)     !Account for A_MLS Symmetry
        A_MLS(i,3,1) = A_MLS(i,1,3)     !Account for A_MLS Symmetry
        A_MLS(i,3,2) = A_MLS(i,2,3)     !Account for A_MLS Symmetry
        do ii = 1,3
          do jj =1,3
            A_input_matrix(ii,jj) = A_MLS(i,ii,jj)
          end do
        end do
        !Amatrix_temp(1,1) = temp  !sum_wab(i)
        i_A_MLS_step = 1
        call LU_decomposition(A_input_matrix,b_input_vector,length,
     &                        A_LU_equivalent_matrix,i_A_MLS_step,
     &                        i_Singular,x_solution,i)      
        !i_Singular = 1
        if(i_Singular.gt.0)then
          !print*,'singular MLS matrix!'
          !print*,'i',i
          !print*,'xp, yp, rhop ',xp(i),yp(i),rhop(i)
          !print*,'sum_wab(i)',sum_wab(i)
          !print*,'Amatrix_temp(1,1)',Amatrix_temp(1,1)
          !stop
          if(A_MLS_11.eq.0.0)then
            write(80,*)'A_MLS_11.eq.0.0'
            stop
          end if
          i_singularMatrix_count = i_singularMatrix_count +1
          beta0_MLS(i) = 1.0/A_MLS_11  !1.0/temp   !1.0/sum_wab(i)
          beta1x_MLS(i)= 0.0
          beta1z_MLS(i)= 0.0
          !print*,'Singular i, xp, zp ',i,xp(i),zp(i)
        else
          !Find inverse matrix by solving system n times 
          !on n successive unit vectors (Gerald & Wheatley, Numerical Analysis) 
          i_A_MLS_step = 2
          do jj = 1,3
            do ii =1,3
              b_input_vector(ii) = 0.0
            end do
            b_input_vector(jj) = 1.0
            call LU_decomposition(A_input_matrix,b_input_vector,length,
     &                            A_LU_equivalent_matrix,i_A_MLS_step,
     &                            i_Singular,x_solution,i)
            do ii =1,3
              A_MLS_inverse(ii,jj) = x_solution(ii)
            end do
          end do
          beta0_MLS(i) = A_MLS_inverse(1,1)
          beta1x_MLS(i)= A_MLS_inverse(2,1)
          beta1z_MLS(i)= A_MLS_inverse(3,1)
                    
        endif
        
      end do
             
c      print*,'itime, MLS i_singularMatrix_count ',
c     &        itime, i_singularMatrix_count
      sum_wab = 0.0
      rhop_sum = 0.0
      i_MLS_part = 2
      
      call ac_MLS   !Perform the corrected density summation
       
      do i=nstart_MLS,np
        rhop_sum_MLS(i) = rhop_sum_MLS(i) + beta0_MLS(i)*adh*pm(i)  !self contribution
        rhop(i) = rhop_sum_MLS(i)
        call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
      enddo


      return
      end

