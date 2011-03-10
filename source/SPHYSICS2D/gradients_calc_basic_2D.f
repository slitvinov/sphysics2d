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


      subroutine gradients_calc(i,j,dux,duz)

      include 'common.2D'  

c     --- Comment out/in the gradients as required ---
      
c     --- Use CSPH terms ---

      Vj = pVol(j)   !=pm(j)/rhop(j)

      frx_norm_i = Vj*frxi
      frz_norm_i = Vj*frzi


      dudx_CSPH(i) = dudx_CSPH(i) - dux*frx_norm_i
      dudz_CSPH(i) = dudz_CSPH(i) - dux*frz_norm_i
      dwdx_CSPH(i) = dwdx_CSPH(i) - duz*frx_norm_i
      dwdz_CSPH(i) = dwdz_CSPH(i) - duz*frz_norm_i
 
      Vi = pVol(i)   !=pm(i)/rhop(i)

      frx_norm_j = Vi*frxj
      frz_norm_j = Vi*frzj

      dudx_CSPH(j) = dudx_CSPH(j) - dux*frx_norm_j
      dudz_CSPH(j) = dudz_CSPH(j) - dux*frz_norm_j
      dwdx_CSPH(j) = dwdx_CSPH(j) - duz*frx_norm_j
      dwdz_CSPH(j) = dwdz_CSPH(j) - duz*frz_norm_j
               
      drho = rhop(i) - rhop(j)
      drhodx_CSPH(i) = drhodx_CSPH(i) - drho*frx_norm_i
      drhodz_CSPH(i) = drhodz_CSPH(i) - drho*frz_norm_i
      drhodx_CSPH(j) = drhodx_CSPH(j) - drho*frx_norm_j
      drhodz_CSPH(j) = drhodz_CSPH(j) - drho*frz_norm_j

      dTE = TEp(i) - TEp(j)
      dTEdx_CSPH(i) = dTEdx_CSPH(i) -dTE*frx_norm_i
      dTEdz_CSPH(i) = dTEdz_CSPH(i) -dTE*frz_norm_i
      dTEdx_CSPH(j) = dTEdx_CSPH(j) -dTE*frx_norm_j
      dTEdz_CSPH(j) = dTEdz_CSPH(j) -dTE*frz_norm_j

      return
      end
