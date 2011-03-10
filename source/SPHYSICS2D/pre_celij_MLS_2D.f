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

      subroutine pre_celij_MLS(j1,j2,kind_p1,ini_kind_p2,lx2)
c
      include 'common.2D'  
      
       do kind_p2=ini_kind_p2,2

        if(nc(j2,kind_p2).ne.0) then

        do ii=1,nc(j1,kind_p1)
          i = ibox(j1,kind_p1,ii)
         
          do jj=1,nc(j2,kind_p2)
           j = ibox(j2,kind_p2,jj)
              
           if(i.ge.nstart_MLS.and.j.ge.nstart_MLS)then

              drx = xp(i) - xp(j)
              drz = zp(i) - zp(j)

              call periodicityCorrection(i,j,drx,drz,lx2)

              rr2 = drx*drx + drz*drz
            
              if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then


                call kernel(drx,drz,i,j,j1,j2,rr2)  
                
                  V_j = pVol(j)  
                  V_i = pVol(i)
                  if(i_MLS_part.eq.1)then   !1st sweep for MLS
                    !--- i-Matrix ---
                    A_MLS(i,1,1) = A_MLS(i,1,1) + Wab*V_j
                    A_MLS(i,1,2) = A_MLS(i,1,2) + drx*Wab*V_j
                    A_MLS(i,1,3) = A_MLS(i,1,3) + drz*Wab*V_j
                    !A_MLS(i,2,1) = A_MLS(i,2,1) + drx*Wab*V_j      !Same as A_MLS(i,1,2)
                    A_MLS(i,2,2) = A_MLS(i,2,2) + drx*drx*Wab*V_j
                    A_MLS(i,2,3) = A_MLS(i,2,3) + drx*drz*Wab*V_j
                    !A_MLS(i,3,1) = A_MLS(i,3,1) + drz*Wab*V_j      !Same as A_MLS(i,1,3)
                    !A_MLS(i,3,2) = A_MLS(i,3,2) + drz*drx*Wab*V_j  !Same as A_MLS(i,2,3)
                    A_MLS(i,3,3) = A_MLS(i,3,3) + drz*drz*Wab*V_j
                    !--- j-Matrix ---
                    A_MLS(j,1,1) = A_MLS(j,1,1) + Wab*V_i
                    A_MLS(j,1,2) = A_MLS(j,1,2) - drx*Wab*V_i
                    A_MLS(j,1,3) = A_MLS(j,1,3) - drz*Wab*V_i
                    !A_MLS(j,2,1) = A_MLS(j,2,1) - drx*Wab*V_i      !Same as A_MLS(j,1,2)
                    A_MLS(j,2,2) = A_MLS(j,2,2) + drx*drx*Wab*V_i
                    A_MLS(j,2,3) = A_MLS(j,2,3) + drx*drz*Wab*V_i
                    !A_MLS(j,3,1) = A_MLS(j,3,1) - drz*Wab*V_i      !Same as A_MLS(j,1,3)
                    !A_MLS(j,3,2) = A_MLS(j,3,2) + drz*drx*Wab*V_i  !Same as A_MLS(j,2,3)
                    A_MLS(j,3,3) = A_MLS(j,3,3) + drz*drz*Wab*V_i
                  else !if(icall.eq.2)then    !2nd Sweep for MLS
                    Wab_MLSi = (beta0_MLS(i) 
     &                        + beta1x_MLS(i)*drx 
     &                        + beta1z_MLS(i)*drz)*Wab
                    Wab_MLSj = (beta0_MLS(j) 
     &                        - beta1x_MLS(j)*drx 
     &                        - beta1z_MLS(j)*drz)*Wab
                    rhop_sum_MLS(i) = rhop_sum_MLS(i) + pm(j)*Wab_MLSi
                    rhop_sum_MLS(j) = rhop_sum_MLS(j) + pm(i)*Wab_MLSj
                  endif   !End of:   if(i_MLS_part.eq.1)then
              
              endif ! if q<2            
             endif   !if(i.ge.nstart_MLS.and.j.ge.nstart_MLS)then 
            enddo ! Box jj
          enddo ! Box ii
        endif  ! Box jj is not empty
      enddo   ! Kind of particle
      
      return
      end
