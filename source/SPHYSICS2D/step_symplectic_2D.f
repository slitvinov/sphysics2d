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


      subroutine step
      include 'common.2D'
      

	
      dt2 = 0.5*dt

c	Predictor- Corrector Method

c  ...  compute acceleration
c

      if(time*1./1000. .eq. int(time*1./1000.) ) then
        write(80,*) 'In step:   time = ',time,' itime = ',itime
        write(*,*) 'In step:   time = ',time,' itime = ',itime
      endif


	call check_limits_2D	!<=========  TEMPORARY
	call ini_divide(2)
	call divide(nbp1,np,2)	


       if(iopt_movingObject.eq.1)then
         call recover_list
         call divide(nbfp1,nb,1)
       endif

      call ac		

 
c
c  ... compute corrections to:
c         (a) rate of change of velocity due to
c             (i)  body forces, and
c             (ii) boundary forces
c
c         (b) rate of change of position (free-surface, XSPH)
c
C
      call correct


      if (ivar_dt.eq.1) then
	    call variable_time_step  
	    ddt_p=dt_new
      endif
	

c
c  ... predictor
c

       do i=nstart,np
         xo(i)   = xp(i)
         zo(i)   = zp(i)
         uo(i)   = up(i)
         wo(i)   = wp(i)
         po(i)   = p(i)
         rhoo(i) = rhop(i)
         TEo(i)=TEp(i)
         
      enddo

      do i=nbp1,np
         rhop(i) = rhoo(i) + dt2*rdot(i)
         pVol(i) = pm(i)/rhop(i)
         TEp(i)=TEo(i)+dt2*TEdot(i)
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
      enddo  


      if(iBC.eq.2)then	!Hughes & Grahams 2009 correction 
         nstep_DBC=nstep_DBC+1
         iDBC=0
         if(mod(nstep_DBC,ndt_DBCPerform).eq.0)then
           do i=1,nb
              rhop(i) = rhoo(i) + dt2*rdot(i)
              if (rhop(i).lt.rho0) rhop(i)=rho0
              pVol(i) = pm(i)/rhop(i)
              TEp(i)=TEo(i)+dt2*TEdot(i)
              call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
           enddo  
           iDBC=1
           nstep_DBC=0
         endif
      endif 

      
      do i=nbp1,np
         xp(i) = xo(i) + dt2*xdot(i)
         zp(i) = zo(i) + dt2*zdot(i)
         
         up(i) = uo(i) + dt2*udot(i)
         wp(i) = wo(i) + dt2*wdot(i)
      
      enddo
           
      
      call movingObjects(0)      

c
c  ...  corrector
c
      call check_limits_2D       !<=========  TEMPORARY
      call ini_divide(2)
      call divide(nbp1,np,2)	

       if(iopt_movingObject.eq.1)then
         call recover_list
         call divide(nbfp1,nb,1)
       endif


      call ac
      call correct


      !-- Perform final integration correction --

      do i=nbp1,np
            epsilon_rdot = -rdot(i)/rhop(i)
            rhop(i)  = rhoo(i) *
     &            (2.0 - epsilon_rdot*dt)/(2.0 + epsilon_rdot*dt)
            pVol(i) = pm(i)/rhop(i)
            !- Thermal Energy -
            TEp(i) = TEp(i) + dt2*TEdot(i)
c  ...      compute pressure using the equation of state
            call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
      enddo

      if(iDBC.eq.1)then	!Hughes & Grahams 2009 correction 
           do i=1,nb
            epsilon_rdot = -rdot(i)/rhop(i)
            rhop(i)  = rhoo(i) *
     &            (2.0 - epsilon_rdot*dt)/(2.0 + epsilon_rdot*dt)
            if (rhop(i).lt.rho0) rhop(i)=rho0
            pVol(i) = pm(i)/rhop(i)
            !- Thermal Energy -
            TEp(i) = TEp(i) + dt2*TEdot(i)
c  ...      compute pressure using the equation of state
            call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
           enddo
      endif 
      
      
      do i=nbp1,np
         up(i)   = uo(i) + dt*udot(i)
         wp(i)   = wo(i) + dt*wdot(i)
         
         xp(i)   = xo(i) + dt2*(uo(i) + up(i))
         zp(i)   = zo(i) + dt2*(wo(i) + wp(i))
         
      enddo

       
      call movingObjects(1)

      if (ivar_dt.eq.1) then
          call variable_time_step 
         ddt_c=dt_new
         dt=min(ddt_p,ddt_c)  !Time step to be used in next loop
      endif
      
      end

