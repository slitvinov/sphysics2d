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

      subroutine movingWedge

      include 'common.2D'
      
      
      t = time

      do i = 2,i_RBend
          
        if(t.gt.RBtime(i-1).and.t.le.RBtime(i))then
          t_factor = (t - RBtime(i-1))/(RBtime(i) - RBtime(i-1))

          do iip = nbfp1,nb
            xp_old = xp(iip)
            zp_old = zp(iip)

            xp(iip) = xp_Wedge_initial(iip-nbf) - RBxpos(i-1)
     &                - t_factor*(RBxpos(i) - RBxpos(i-1))

            zp(iip) = zp_Wedge_initial(iip-nbf)
     &                -(xp_Wedge_initial(iip-nbf) - xp(iip))*bslope
            up(iip) = (xp(iip) - xp_old)/dt
            wp(iip) = (zp(iip) - zp_old)/dt
            

          end do
        end if
      
      end do

	      
	return      
      end