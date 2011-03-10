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

c
c	It calculates the Energy (Thermal, Kin & Pot)
C	Both for boundary particles and fluid particles
C

      subroutine energy(g)
      include 'common.2D'
      
c


      Eki_b=0. ! Boundary particles
      Epo_b=0.
      TE_b=0.

      do i=nstart,nb
           Eki_b =Eki_b+0.5*pm(i)*(up(i)*up(i)
     +		                    +wp(i)*wp(i))
           Epo_b=Epo_b+g*pm(i)*zp(i)
           TE_b=TE_b+TEp(i)*pm(i)
      enddo

      Eki_p=0.        ! Fluid Particles
      Epo_p=0.
      TE_p=0.

      do i=nbp1,np
           Eki_p =Eki_p+0.5*pm(i)*(up(i)*up(i)
     +		                    +wp(i)*wp(i))
           Epo_p=Epo_p+g*pm(i)*zp(i)
           TE_p=TE_p+TEp(i)*pm(i)
      enddo

      write(50,100) time,Eki_p,Epo_p,TE_p,Eki_b,Epo_b,TE_b

100   format(7e16.8)


      return 
      end
