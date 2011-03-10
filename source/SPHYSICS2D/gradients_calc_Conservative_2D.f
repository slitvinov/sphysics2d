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


      subroutine gradients_calc(i,j,drho,drhou,drhow,dTotalE)

      include 'common.2D'  
      
c     --- Comment out/in the gradients as required ---
      
c     --- Use CSPH terms ---

      Vj = pVol(j)   !=pm(j)/rhop(j)

      frx_norm_i = Vj*frxi
      frz_norm_i = Vj*frzi

c      dudx_CSPH(i) = dudx_CSPH(i) - dux*frx_norm_i
cc     dudy_CSPH(i) =
c      dudz_CSPH(i) = dudz_CSPH(i) - dux*frz_norm_i
cc     dvdx_CSPH(i) =
cc     dvdy_CSPH(i) =
cc     dvdz_CSPH(i) =
c      dwdx_CSPH(i) = dwdx_CSPH(i) - duz*frx_norm_i
cc     dwdy_CSPH(i) =
c      dwdz_CSPH(i) = dwdz_CSPH(i) - duz*frz_norm_i
 
      Vi = pVol(i)   !=pm(i)/rhop(i)

      frx_norm_j = Vi*frxj
      frz_norm_j = Vi*frzj

c      dudx_CSPH(j) = dudx_CSPH(j) - dux*frx_norm_j
cc     dudy_CSPH(i) =
c      dudz_CSPH(j) = dudz_CSPH(j) - dux*frz_norm_j
cc     dvdx_CSPH(i) =
cc     dvdy_CSPH(i) =
cc     dvdz_CSPH(i) =
c      dwdx_CSPH(j) = dwdx_CSPH(j) - duz*frx_norm_j
cc     dwdy_CSPH(i) =
c      dwdz_CSPH(j) = dwdz_CSPH(j) - duz*frz_norm_j
               
      drhodx_CSPH(i) = drhodx_CSPH(i) - drho*frx_norm_i
c     drhody_CSPH(i) = drhody_CSPH(i) - drho*fry_norm_i
      drhodz_CSPH(i) = drhodz_CSPH(i) - drho*frz_norm_i

      drhoudx_CSPH(i) = drhoudx_CSPH(i) - drhou*frx_norm_i
c     drhoudy_CSPH(i) = drhoudy_CSPH(i) - drhou*fry_norm_i
      drhoudz_CSPH(i) = drhoudz_CSPH(i) - drhou*frz_norm_i
c     drhovdx_CSPH(i) = drhovdx_CSPH(i) - drhov*frx_norm_i
c     drhovdy_CSPH(i) = drhovdy_CSPH(i) - drhov*fry_norm_i
c     drhovdz_CSPH(i) = drhovdz_CSPH(i) - drhov*frz_norm_i
      drhowdx_CSPH(i) = drhowdx_CSPH(i) - drhow*frx_norm_i
c     drhowdy_CSPH(i) = drhowdy_CSPH(i) - drhow*fry_norm_i
      drhowdz_CSPH(i) = drhowdz_CSPH(i) - drhow*frz_norm_i
               
      drhodx_CSPH(j) = drhodx_CSPH(j) - drho*frx_norm_j
c     drhody_CSPH(j) = drhody_CSPH(j) - drho*fry_norm_j
      drhodz_CSPH(j) = drhodz_CSPH(j) - drho*frz_norm_j
      drhoudx_CSPH(j) = drhoudx_CSPH(j) - drhou*frx_norm_j
c     drhoudy_CSPH(j) = drhoudy_CSPH(j) - drhou*fry_norm_j
      drhoudz_CSPH(j) = drhoudz_CSPH(j) - drhou*frz_norm_j
c     drhovdx_CSPH(j) = drhovdx_CSPH(j) - drhov*frx_norm_j
c     drhovdy_CSPH(j) = drhovdy_CSPH(j) - drhov*fry_norm_j
c     drhovdz_CSPH(j) = drhovdz_CSPH(j) - drhov*frz_norm_j
      drhowdx_CSPH(j) = drhowdx_CSPH(j) - drhow*frx_norm_j
c     drhowdy_CSPH(j) = drhowdy_CSPH(j) - drhow*fry_norm_j
      drhowdz_CSPH(j) = drhowdz_CSPH(j) - drhow*frz_norm_j

      dTotalEdx_CSPH(i) = dTotalEdx_CSPH(i) -dTotalE*frx_norm_i
c     dTotalEdy_CSPH(i) = dTotalEdy_CSPH(i) -dTotalE*fry_norm_i
      dTotalEdz_CSPH(i) = dTotalEdz_CSPH(i) -dTotalE*frz_norm_i
      dTotalEdx_CSPH(j) = dTotalEdx_CSPH(j) -dTotalE*frx_norm_j
c     dTotalEdy_CSPH(j) = dTotalEdy_CSPH(j) -dTotalE*fry_norm_j
      dTotalEdz_CSPH(j) = dTotalEdz_CSPH(j) -dTotalE*frz_norm_j



      return
      end
