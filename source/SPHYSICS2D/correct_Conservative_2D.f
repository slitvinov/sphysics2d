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

      subroutine correct
c
      include 'common.2D'
      

c
c     ...  account for body forces
c
      do i=nbp1,np

        !udot(i) = udot(i) + grx*iflag(i)
        !wdot(i) = wdot(i) + grz*iflag(i)


        aVol(i)    = pVol(i)*aVol(i)
        ar(i)      = pVol(i)*ar(i)
        ax(i)      = pVol(i)*ax(i) + Volrho(i)*grx*iflag(i)
        az(i)      = pVol(i)*az(i) + Volrho(i)*grz*iflag(i)

        aTE(i)  = pVol(i)*aTE(i) 

      enddo

c
c      ...  account for XSPH
c

      do i=nbp1,np

        xdot(i) = up(i) + xcor(i)
        zdot(i) = wp(i) + zcor(i)

        
      enddo


      return
      end
