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

c            --- Viscous Diffusion terms (Morris et al 1997) ---    
 
      subroutine viscosity(dot,drx,drz,dux,duz,rr2,
     +           cbar,robar,one_over_rhobar,i,j,j1,j2,term2i,term2j)

      include 'common.2D'  
      

c         --- Viscous Diffusion terms (Lo & Shao 2002) ---

      tempi =2.0*viscos_val*one_over_rhobar*
     &     ((drx*frxi+drz*frzi)/(rr2 + eta2))

      tempj =2.0*viscos_val*one_over_rhobar*
     &     ((drx*frxj+drz*frzj)/(rr2 + eta2))

      ax(i) = ax(i) + pm(j)*tempi*dux
      az(i) = az(i) + pm(j)*tempi*duz

      ax(j) = ax(j) - pm(i)*tempj*dux
      az(j) = az(j) - pm(i)*tempj*duz

      term2i=  -0.5 * tempi *( dux*dux+duz*duz)
      term2j=  -0.5 * tempj *( dux*dux+duz*duz)

      return
      end
