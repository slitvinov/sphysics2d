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
	
      dt=min(ddt_p,ddt_c) 
      dt2 = 0.5*dt
      dt_sq = dt*dt

c	  BEEMAN's SOLVER

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
C     STORE VALUES AT t-dt,t
      do i=nstart,np
         udotm1(i)=udoto(i)
         wdotm1(i)=wdoto(i)
         rdotm1(i)=rdoto(i)
         xo(i)   = xp(i)
         zo(i)   = zp(i)
         uo(i)   = up(i)
         wo(i)   = wp(i)
         po(i)   = p(i)
         rhoo(i) = rhop(i)
         TEo(i)=TEp(i)
      enddo
      
      call ac		

      call correct

      if (ivar_dt.eq.1) then
        call variable_time_step  
        ddt_p=dt_new
      endif

c
c  ... predictor
c     THE VARIABLES WITH SUFFIXES o STORE VALUES AT TIME t
      
      do i=nbp1,np
         xp(i) = xo(i) +
     &           xdot(i)*dt+2.*udot(i)*dt_sq/3.-udotm1(i)*dt_sq/6.
         zp(i) = zo(i) + 
     &           zdot(i)*dt+2.*wdot(i)*dt_sq/3.-wdotm1(i)*dt_sq/6.
 
         up(i) = uo(i) + 1.5*udot(i)*dt-0.5*udotm1(i)*dt
         wp(i) = wo(i) + 1.5*wdot(i)*dt-0.5*wdotm1(i)*dt

         udoto(i)=udot(i)
         wdoto(i)=wdot(i)
      enddo

      do i=nbp1,np
         rhop(i) = rhoo(i) + 1.5*rdot(i)*dt-0.5*rdotm1(i)*dt   
         pVol(i) = pm(i)/rhop(i)  
c  ...   compute pressure using the equation of state
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
c  ...   Compute thermal energy
         TEp(i)=TEo(i)+dt2*TEdot(i)
         rdoto(i)=rdot(i)
      enddo
 

      if(iBC.eq.2)then	!Hughes & Grahams 2009 correction 
         nstep_DBC=nstep_DBC+1
         iDBC=0
         if(mod(nstep_DBC,ndt_DBCPerform).eq.0)then
           do i=1,nb
             rhop(i) = rhoo(i) + 1.5*rdot(i)*dt-0.5*rdotm1(i)*dt   
             if (rhop(i).lt.rho0) rhop(i)=rho0
             pVol(i) = pm(i)/rhop(i)  
             call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
             TEp(i)=TEo(i)+dt2*TEdot(i)
             rdoto(i)=rdot(i)
           enddo  
           iDBC=1
           nstep_DBC=0
         endif
      endif  
   
     
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

      if (ivar_dt.eq.1) then
	    call variable_time_step 
	    ddt_c=dt_new
      endif

      do i=nbp1,np
         xp(i) = xo(i) + 
     +           dt*xdot(i)+udot(i)*dt_sq/6.+udoto(i)*dt_sq/3.
         zp(i) = zo(i) +
     +           dt*zdot(i)+wdot(i)*dt_sq/6.0+wdoto(i)*dt_sq/3.
         up(i) = uo(i) + 
     +           (5.0*udot(i)+8.0*udoto(i)-udotm1(i))*dt/12.
         wp(i) = wo(i) + 
     +           (5.0*wdot(i)+8.0*wdoto(i)-wdotm1(i))*dt/12.   
      enddo


      do i=nbp1,np
         rhop(i) = rhoo(i) +
     +            (5.0*rdot(i)+8.0*rdoto(i)-rdotm1(i))*dt/12.0
         pVol(i) = pm(i)/rhop(i)    
c  ...   compute pressure using the equation of state
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
c  ...   Compute thermal energy
         TEp(i)=TEo(i)+dt2*TEdot(i)  
      enddo

      if(iDBC.eq.1)then	!Hughes & Grahams 2009 correction 
           do i=1,nb
            rhop(i) = rhoo(i) +
     +            (5.0*rdot(i)+8.0*rdoto(i)-rdotm1(i))*dt/12.0
            if (rhop(i).lt.rho0) rhop(i)=rho0
            pVol(i) = pm(i)/rhop(i)    
            call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
            TEp(i)=TEo(i)+dt2*TEdot(i)
           enddo
      endif 


      call movingObjects(1)

      end
c