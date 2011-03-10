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
c*********************************************************************
c
      subroutine movingPaddle(i_half_dt)

      include 'common.2D'
      
     
      double precision xb_min_double, xb_max_double, xl_double
      double precision zb_min_double, zb_max_double, zl_double
      
c     i_half_dt = 1,2  depending on whether subroutine is being called 
c                      during the first/second half of the timestep


      if(i_paddleType.eq.1)then
c       ---  Piston paddle ---                      
           
        if (time.gt.(twinitial(ind_wm)+twinitial(ind_wm+1))*0.5) then

          ind_wm=ind_wm+1
          ind1_wm=ind_wm+1
        endif


c       twopi_o_Period=twopi/Period(ind_wm)
c       twopi_o_Period1=twopi/Period(ind1_wm)

c       write(*,*), time, twopi_o_Period, twopi_o_Period1

c       wm_coeff1=(1-tanh(time-twinitial(ind_wm)))
c       wm_coeff2=(1+tanh(time-twinitial(ind_wm)))

       if(Nfreq.gt.1)then

        do i=nwavemaker_ini,nwavemaker_end

          xp(i)=xp_ini(i)+
     +     0.5*(A_wavemaker(ind_wm)*(1-tanh(time-twinitial(ind_wm)))*
     +        sin(2*pi*(time-twinitial(ind_wm))/Period(ind_wm))+
     +        A_wavemaker(ind1_wm)*(1+tanh(time-twinitial(ind_wm)))*
     +        sin(2*pi*(time-twinitial(ind1_wm))/Period(ind1_wm)))

          up(i)=0.5*2*pi*((A_wavemaker(ind_wm)/Period(ind_wm))*
     +        (1-tanh(time-twinitial(ind_wm)))*
     +        cos(2*pi*(time-twinitial(ind_wm))/Period(ind_wm))+
     +        (A_wavemaker(ind1_wm)/Period(ind1_wm))*
     +        (1+tanh(time-twinitial(ind_wm)))*
     +        cos(2*pi*(time-twinitial(ind1_wm))/Period(ind1_wm)))

        enddo
       else
         wave_freq = 2.0*pi/Period(1)
        do i=nwavemaker_ini,nwavemaker_end

          xp(i)=xp_ini(i)+
     +        A_wavemaker(1)*sin(wave_freq*time)

          up(i)=wave_freq*A_wavemaker(1)*cos(wave_freq*time)

        enddo
         
       endif

      else if(i_paddleType.eq.2)then
c         ---  Piston-flap paddle ---             
        t = time          
        t_frac = 1.0   !t/Tr
        t_factor = 1.0
        !- Optional Ramp-up Function -
        !if(t_frac.lt.1.0)then
        !  t_factor = (-2.0*(t_frac**3) + 3.0*(t_frac**2))
        !else
        !  t_factor = 1.0
        !end if
        wave_freq = 2.0*pi/Period(1)
        
        do i = nwavemaker_ini,nwavemaker_end
          z_factor = 1.0 + (zp(i)-paddle_SWL)/
     &               (paddle_SWL + flap_length)
          xp(i) = X_PaddleCentre
     &            + z_factor*t_factor*0.5*stroke*cos(wave_freq*t)
          up(i) = -wave_freq
     &             *z_factor*t_factor*0.5*stroke*sin(wave_freq*t)   
        end do

c       yp(iip) remains unchanged
c       zp(iip) remains unchanged
c       vp(iip) remains unchanged
c       wp(iip) remains unchanged

      else if(i_paddleType.eq.3)then
c         ---  Piston paddle with Prescribed Motion ---             
        t = time          
        do i = 2,n_paddleData
          if(t.gt.time_paddle(i-1).and.t.le.time_paddle(i))then
            t_factor = (t - time_paddle(i-1))/
     &               (time_paddle(i) - time_paddle(i-1))
        
            do ii = nwavemaker_ini,nwavemaker_end
              xp(ii) = x_paddle(i-1)
     &               + t_factor*(x_paddle(i)-x_paddle(i-1))
              up(ii) = u_paddle(i-1)
     &               + t_factor*(u_paddle(i)-u_paddle(i-1))
            enddo
          endif   
        end do
      end if
c       zp(iip) remains unchanged
c       wp(iip) remains unchanged


c     Normals are now calculated in subroutine updateNormals
      if(iBC.eq.1)then   !Repulsive Force BC
        call updateNormals_2D(nbfp1,nb)
      end if
      
	return      
      end
c
c*********************************************************************
c
