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

      subroutine pre_celij_Shepard(j1,j2,kind_p1,ini_kind_p2,lx2)
c
      include 'common.2D'  
      
       do kind_p2=ini_kind_p2,2

        if(nc(j2,kind_p2).ne.0) then

        do ii=1,nc(j1,kind_p1)
          i = ibox(j1,kind_p1,ii)
         
          do jj=1,nc(j2,kind_p2)
           j = ibox(j2,kind_p2,jj)

            if(i.ge.nstart_Shepard.and.j.ge.nstart_Shepard)then      
		          
              drx = xp(i) - xp(j)
              drz = zp(i) - zp(j)

              call periodicityCorrection(i,j,drx,drz,lx2)

              rr2 = drx*drx + drz*drz
            
              if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then

                call kernel(drx,drz,i,j,j1,j2,rr2)  
                 
                !- Debug code -
                i_particleCheck = -23197
                i_rankCheck = 0
                dpN = 0.01*vnorm_mass
599             format(a8,2(1x,i7),i5,2(1x,i2),f12.6)             
                   !Place in separate ifs            
                  !Place in separate ifs            

c     Sum_wab to calculate the normalized kernel
			if(j.ge.nstart_Shepard)then
			  V_j = pVol(j)  !=pm(j)/rhop(j)               
			  sum_wab(i)  = sum_wab(i)  + Wab*V_j   
			  rhop_sum(i) = rhop_sum(i) + pm(j)*Wab       
                 if(0.eq.i_rankCheck.and.i.eq.i_particleCheck)then
                    idrx = nint(drx/dpN)
                    idrz = nint(drz/dpN)
                    write(*,599)
     &              'pC:I,j ',i,j,lx2,idrx,idrz,Wab*V_j
                  endif 
			end if

			if(i.ge.nstart_Shepard)then                   
			  V_i = pVol(i)   !=pm(i)/rhop(i)
			  sum_wab(j)  = sum_wab(j)  + Wab*V_i 
			  rhop_sum(j) = rhop_sum(j) + pm(i)*Wab       
                 if(0.eq.i_rankCheck.and.j.eq.i_particleCheck)then
                    idrx = -nint((drx)/dpN)
                    idrz = -nint((drz)/dpN)
                    write(*,599)
     &              'pC:i,J ',i,j,lx2,idrx,idrz,Wab*V_i
                  endif 
			endif
            
              endif ! if q<2 
             endif   !if(i.ge.nstart_Shepard.and.j.ge.nstart_Shepard)           
          enddo ! Box jj
        enddo ! Box ii
        endif  ! Box jj is not empty
      enddo   ! Kind of particle
      
      return
      end
