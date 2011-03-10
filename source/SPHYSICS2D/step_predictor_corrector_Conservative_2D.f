c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Alejandro Crespo, Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
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

         Volo(i)       = pVol(i)
         Volrhoo(i)    = Volrho(i)
         Volrhouo(i)   = Volrhou(i)
         Volrhowo(i)   = Volrhow(i)
         VolTotalEo(i)  = VolTotalE(i)
         
      enddo

      do i=nstart,np
         !rhop(i) = rhoo(i) + dt2*rdot(i)
         !TEp(i)=TEo(i)+dt2*TEdot(i)
c  ...   Compute volume and density
         pVol(i) = Volo(i) + dt2*aVol(i)
         Volrho(i) = Volrhoo(i) + dt2*ar(i)
         rhop(i) = Volrho(i) / pVol(i)
      end do
      
      do i=nbp1,np
         xp(i) = xo(i) + dt2*xdot(i)
         zp(i) = zo(i) + dt2*zdot(i)

         Volrhou(i) = Volrhouo(i) + dt2*ax(i)
         Volrhow(i) = Volrhowo(i) + dt2*az(i)
         up(i) = Volrhou(i) /Volrho(i)
         wp(i) = Volrhow(i) /Volrho(i)
      enddo
      
      do i=nstart,np
c  ...   Compute Total and thermal energy
         VolTotalE(i) = VolTotalEo(i) + dt2*aTE(i)
         TotalE(i) = VolTotalE(i) / pVol(i)
         TEp(i) = (TotalE(i)/rhop(i)) - 0.5*(up(i)*up(i) + wp(i)*wp(i))
c  ...   compute pressure using the equation of state
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
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



      if (ivar_dt.eq.1) then
	    call variable_time_step 
	    ddt_c=dt_new
           dt=min(ddt_p,ddt_c)  !Time step to be used in next loop
      endif

      do i=nstart,np
         !rhop(i) = rhoo(i) + dt2*rdot(i)
         !TEp(i)=TEo(i)+dt2*TEdot(i)
         pVol(i) = Volo(i) + dt2*aVol(i)
         Volrho(i) = Volrhoo(i) + dt2*ar(i)
c  ...   Compute Total energy
         VolTotalE(i) = VolTotalEo(i) + dt2*aTE(i)
      end do

      do i=nbp1,np
         xp(i) = xo(i) + dt2*xdot(i)
         zp(i) = zo(i) + dt2*zdot(i)

         Volrhou(i) = Volrhouo(i) + dt2*ax(i)
         Volrhow(i) = Volrhowo(i) + dt2*az(i)
      enddo

      
      !-- Perform final integration correction --
      do i=nstart,np
         !rhop(i) = 2.*rhop(i) - rhoo(i)
         !TEp(i)=2.*TEp(i)-TEo(i)
         pVol(i)   = 2.*pVol(i) - Volo(i)
         Volrho(i) = 2.*Volrho(i) - Volrhoo(i)
         rhop(i) = Volrho(i) / pVol(i)
      end do
      
      do i=nbp1,np
         
         xp(i)   = 2.*xp(i) - xo(i)
         zp(i)   = 2.*zp(i) - zo(i)
         
         Volrhou(i)   = 2.*Volrhou(i) - Volrhouo(i)
         Volrhow(i)   = 2.*Volrhow(i) - Volrhowo(i)
         up(i) = Volrhou(i) / Volrho(i)
         wp(i) = Volrhow(i) / Volrho(i)
      enddo

      do i=nstart,np
c  ...   Compute thermal energy
         VolTotalE(i)=2.*VolTotalE(i)-VolTotalEo(i)
         TotalE(i) = VolTotalE(i) / pVol(i)
         TEp(i) = TotalE(i)/rhop(i) - 0.5*(up(i)*up(i) + wp(i)*wp(i))
c  ...   compute pressure using the equation of state
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
      enddo
       
      call movingObjects(1)

      end
