c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr Alejandro Crespo, Dr Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "

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


      subroutine vorticity_calc(i,j,dux,duz)

      include 'common.2D'  	
      
      pmj_over_rhoi = pm(j)/rhop(i)

      curl3y_a(i)= curl3y_a(i)+duz*frxi*pmj_over_rhoi
      curl3y_b(i)= curl3y_b(i)+dux*frzi*pmj_over_rhoi

      pmi_over_rhoj = pm(i)/rhop(j)

      curl3y_a(j)= curl3y_a(j)+duz*frxj*pmi_over_rhoj
      curl3y_b(j)= curl3y_b(j)+dux*frzj*pmi_over_rhoj


      return
      end
