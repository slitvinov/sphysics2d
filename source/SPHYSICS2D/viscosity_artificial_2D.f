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


c  ...  viscous term from (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.1b)
c
      
      subroutine viscosity(dot,drx,drz,dux,duz,rr2,
     +           cbar,robar,one_over_rhobar,i,j,j1,j2,term2i,term2j)

      include 'common.2D'  	
      


      if (dot.le.0.)then         ! aproaching
         amubar = h*dot/(rr2 + eta2)
         pi_visc = - viscos_val * cbar * amubar / robar
      else                                     ! going away
         pi_visc = 0.
      endif
 
      ax(i) = ax(i) - pm(j) * pi_visc * frxi
      az(i) = az(i) - pm(j) * pi_visc * frzi

      ax(j) = ax(j) + pm(i) * pi_visc * frxj
      az(j) = az(j) + pm(i) * pi_visc * frzj

      term2i=0.5 * pi_visc *( frxi*dux+frzi*duz)
      term2j=0.5 * pi_visc *( frxj*dux+frzj*duz)
	
      return
      end
