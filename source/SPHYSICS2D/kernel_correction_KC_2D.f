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

      subroutine kernel_correction(i,j)

      include 'common.2D'   

c     Sum_wab to calculate the normalized kernel
	if(j.gt.nb)then     
	  V_j = pVol(j)  !=pm(j)/rhop(j)               
	  sum_wab(i)  = sum_wab(i)  + Wab*V_j   
	  rhop_sum(i) = rhop_sum(i) + pm(j)*Wab       
	end if

	if(i.gt.nb)then        
	  V_i = pVol(i)   !=pm(i)/rhop(i)
	  sum_wab(j)  = sum_wab(j)  + Wab*V_i 
	  rhop_sum(j) = rhop_sum(j) + pm(i)*Wab       
	endif

      frxi=frx
      frzi=frz

      frxj=frx
      frzj=frz


      return 
      end
