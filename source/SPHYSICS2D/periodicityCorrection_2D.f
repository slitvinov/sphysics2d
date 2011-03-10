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


      subroutine periodicityCorrection(i,j,drx,drz,lx2)

      include 'common.2D'

      double precision dxp_dble1,dxp_dble2

      !- X-Periodicity Check -
      if(lx2.gt.ncx)then   !if(i_periodicOBs(1).gt.0)then
        if(drx.gt.four_h)then  !if(dry.gt.4.0*h)then
          dxp_dble1 = xmax_container_double-dble(xp(i))
          dxp_dble2 = dble(xp(j))-xmin_container_double
          drx = -real(dxp_dble1+dxp_dble2)
        endif
      elseif(lx2.lt.3)then
        if(drx.lt.-four_h)then  !if(dry.gt.4.0*h)then
          dxp_dble1 = dble(xp(i)) - xmin_container_double
          dxp_dble2 = xmax_container_double - dble(xp(j))
          drx =  real(dxp_dble1+dxp_dble2)
        end if
      end if

      return
      end