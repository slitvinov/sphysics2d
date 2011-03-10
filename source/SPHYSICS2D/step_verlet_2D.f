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
      

C	Verlet Method

c  ...  compute acceleration
c
C
      if(time*1./1000 .eq. int(time*1./1000) ) then
         write(*,*) 'In step:   time = ',time
      endif


      call check_limits_2D
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


C     To calculate time step 
      if (ivar_dt.eq.1) then
         call variable_time_step	  
         dt=dt_new
      endif	
      dtsq_05 = 0.5*dt*dt
      two_dt=2*dt

      do i=nstart,np                                               

          um1(i)  = uo(i)

          wm1(i)  = wo(i)
          xo(i)   = xp(i)

          zo(i)   = zp(i)
          uo(i)   = up(i)

          wo(i)   = wp(i)
          rhom1(i)= rhoo(i)
          rhoo(i) = rhop(i)
          TEm1(i)  = TEo(i)    
          TEo(i) = TEp(i)
      enddo


c     determine which velocity term to use


      nstep_verlet=nstep_verlet+1
      imix=0
      if(mod(nstep_verlet,ndt_VerletPerform).eq.0)then
        imix=1     
        nstep_verlet = 0
      end if
         
      !print*,'itime, nstep_verlet, imix ',itime, nstep_verlet, imix

      if (imix.eq.0) then
        do i=nbp1,np                          
           xp(i) = xo(i) +xdot(i)*dt +udot(i)*dtsq_05        
           zp(i) = zo(i) +zdot(i)*dt +wdot(i)*dtsq_05 
           
           up(i) = um1(i) + two_dt*udot(i)
           wp(i) = wm1(i) + two_dt*wdot(i)
        enddo
      else
        do i=nbp1,np 
           xp(i) = xo(i) +xdot(i)*dt +udot(i)*dtsq_05       
           zp(i) = zo(i) +zdot(i)*dt +wdot(i)*dtsq_05 

           up(i)=uo(i)+udot(i)*dt
           wp(i)=wo(i)+wdot(i)*dt
        enddo
      endif

      if(imix.eq.0) then
        do i=nstart,np                                
           rhop(i) = rhom1(i) + two_dt*rdot(i)
           pVol(i) = pm(i)/rhop(i)
           TEp(i) = TEm1(i) + two_dt*TEdot(i) 
c  ...     compute pressure using the equation of state
           !p(i) = B * ( (rhop(i)/rho0)**i_gamma - 1.)
           call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
        enddo
      else
        do i=nstart,np     
           rhop(i) = rhoo(i) + dt*rdot(i)
           pVol(i) = pm(i)/rhop(i)
           TEp(i) = TEo(i) + dt*TEdot(i) 
c  ...     compute pressure using the equation of state
           !p(i) = B * ( (rhop(i)/rho0)**i_gamma - 1.)
           call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
        enddo
      endif


      call movingObjects(1)

      end