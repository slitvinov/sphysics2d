c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr. Alejandro Crespo, Dr. Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
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

      subroutine movingObjects(i_midPoint)
c
      include 'common.2D'
      
      double precision xb_min_double, xb_max_double, xl_double
      double precision zb_min_double, zb_max_double, zl_double


      if(i_midPoint.gt.0)then     
	
        if(iopt_gate.eq.1)  call movingGate	 
        if(iopt_wavemaker.eq.1)      call movingPaddle(0) !arguement is an integer
        if(iopt_RaichlenWedge.eq.1)  call movingWedge
      
      endif
      
      if(iopt_FloatingBodies.eq.1.and.itime.gt.0)then
        call rigid_body_motion(i_midPoint)
      endif

   
      return
      end

