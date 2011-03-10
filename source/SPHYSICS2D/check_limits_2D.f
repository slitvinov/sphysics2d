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


      subroutine check_limits_2D
      include 'common.2D'
      
      double precision dxp_dble

      zmax=0
      ncases=0

      do i=nbp1,np ! Check if any fluid particle is over the box
        if(zp(i).gt.zmax) zmax=zp(i)
          if(iflag(i).eq.1.and.(xp(i).gt.xmax_container.or.
     +      xp(i).lt.xmin_container.or.zp(i).lt.zmin_container)) then	
            if(i_periodicOBs(1).eq.1.and.i.gt.nb
     +              .and.zp(i).gt.zmin_container)then  !Replace Particle on other side of domain
              if(xp(i).gt.xmax_container)then
                dxp_dble = dble(xp(i)) - xmax_container_double 
                xp(i) = real(dxp_dble +xmin_container_double)
              else if(xp(i).lt.xmin_container)then
                dxp_dble = xmin_container_double - dble(xp(i))   
                xp(i) = real(xmax_container_double - dxp_dble)
              end if
            else
            write(80,*) 'Particle out of limits: ',i, 
     +      '  last X position: ',xp(i), '  last Z-position: ',zp(i)
            write(*,*) 'Particle out of limits: ',i, 
     +      '  last X position: ',xp(i), '  last Z-position: ',zp(i)
             iflag(i)=0
             ncases=ncases+1
             xp(i)=xtrash
             zp(i)=ztrash
             up(i)=0.
             wp(i)=0.
             um1(i)=0.
             uo(i)=0.
             wm1(i)=0.
             wo(i)=0.
            endif
         endif   
	enddo

       !Check to see if Moving BPs moves across periodic boundaries
       if(i_periodicOBs(1).eq.1)then
         do i=nbf+1,nb
           if(zp(i).gt.zmax) zmax=zp(i)
           if(xp(i).gt.xmax_container.or.xp(i).lt.xmin_container)then
             if(xp(i).gt.xmax_container)then
               dxp_dble = dble(xp(i)) - xmax_container_double 
               xp(i) = real(dxp_dble +xmin_container_double)
             else if(xp(i).lt.xmin_container)then
               dxp_dble = xmin_container_double - dble(xp(i))   
               xp(i) = real(xmax_container_double - dxp_dble)
             end if
           endif
         enddo
       endif

      differ=(zmax-zcontrol)
   
       if (differ.gt.0) then     
          nn=int(differ*one_over_2h)+1
          zmax = zcontrol + two_h*nn
          ncz = int( (zmax-zmin)*one_over_2h ) + 1
          nct = ncx*ncz !2D
              do i=nct_ini+1,nct
                 nc(i,1)  = 0 !No Boundary Particles in these cells               
              enddo
       else
         zmax=zcontrol
          nct=nct_ini
          ncz=ncz_ini
       endif  


      if (ncases.ne.0) write(80,*) 'Number ',ncases,' time',time
      if (ncases.ne.0) write(*,*) 'Number ',ncases,' time',time

      if (ncases.ge.100) stop

       if(nct.gt.nct_max.or.ncz.lt.1)then
         write(80,*) ' '
         write(80,*)'ERROR in check_limits_2D.f'
         write(80,*)'nct.gt.nct_max.or.ncz.lt.1'
         write(80,*)'ncx/z    = ',ncx,ncz
         write(80,*)'nct= ',nct
         write(80,*)'nct_max ',nct_max
         write(80,*)'itime  ', itime
c        -- screen printout ---              
         write(*,*) ' '
         write(*,*)'ERROR in check_limits_2D.f'
         write(*,*)'nct.gt.nct_max.or.ncz.lt.1'
         write(*,*)'ncx/z    = ',ncx,ncz
         write(*,*)'nct= ',nct
         write(*,*)'nct_max ',nct_max
         write(*,*)'itime  ', itime
         stop
       endif

      return
	end
