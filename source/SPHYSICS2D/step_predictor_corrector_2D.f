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
      
	
       character supp*4,name*40,detsupp*4
      
      dt2 = 0.5*dt

c	Predictor- Corrector Method

c  ...  compute acceleration
c
      !-- Debug Code------
      itime_check = 10e8   !0    !10e8   !0    !
      i_rank_check = 0
      i_ParticleCheck = 23197  !2366  !nb_local + 1
      

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

      !if(itime.eq.0) call densityFilter
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

c      i = 1087	
c      print*,'After CORRECT, predictor, itime = ',itime
c      print*,'xp(i), zp(i) ',xp(i), zp(i)
c      print*,'ax(i), az(i) ',ax(i), az(i)
c      print*,'dt, dt2 ',dt, dt2 
c      i = 565	
c      print*,'xp(i), zp(i) ',xp(i), zp(i)
c      print*,'ax(i), az(i) ',ax(i), az(i)
c      print*,'dt, dt2 ',dt, dt2 
c      if(itime.eq.3)stop

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
      end do


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
      
      !-- Debug Code------
      if(itime.ge.itime_check.and.0.eq.i_rank_check)then
        print*
        print*,'AFTER CORRECT, PREDICTOR, before movingOBJECTS'
        print*,'rank_id, itime ',rank_id, itime
        i = i_ParticleCheck  
        print*,'i = ',i
        print*,'xp, zp  ',xp(i),zp(i)
        print*,'up, wp  ',up(i),wp(i)
        print*,'rhop, p,  pVol, pm ',rhop(i),p(i) ,pVol(i),pm(i)
        print*,'rhoo, po, Volo, pm ',rhoo(i),po(i),Volo(i),pm(i)
        print*,'TEp ',TEp(i)
        print*,'ar(i) ',ar(i)
        print*,'ax(i) ',ax(i)
        print*,'az(i) ',az(i)
        print*,'pVol(i), Volo(i) ',pVol(i), Volo(i)
        print*,'dudx_CSPH(i)',dudx_CSPH(i)
        print*,'dudz_CSPH(i)',dudz_CSPH(i)
c
        print*,'dwdx_CSPH(i)',dwdx_CSPH(i)
        print*,'dwdz_CSPH(i)',dwdz_CSPH(i)
        print*,'drhodx_CSPH(i)',drhodx_CSPH(i)
        print*,'drhodz_CSPH(i)',drhodz_CSPH(i)
c
        print*,'sum_wab(i) ',sum_wab(i)
       !stop
      end if
      !---------------------

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

      do i=nbp1,np
         rhop(i) = rhoo(i) + dt2*rdot(i)
         TEp(i)=TEo(i)+dt2*TEdot(i)
      end do


      if(iDBC.eq.1)then	!Hughes & Grahams 2009 correction 
           do i=1,nb
              rhop(i) = rhoo(i) + dt2*rdot(i)
              if (rhop(i).lt.rho0) rhop(i)=rho0
              TEp(i)=TEo(i)+dt2*TEdot(i)
           enddo
      endif 

      do i=nbp1,np
         xp(i) = xo(i) + dt2*xdot(i)
         zp(i) = zo(i) + dt2*zdot(i)
         
         up(i) = uo(i) + dt2*udot(i) 
         wp(i) = wo(i) + dt2*wdot(i)

      enddo

      
      !-- Perform final integration correction --
      do i=nbp1,np
         rhop(i) = 2.*rhop(i) - rhoo(i)
         pVol(i) = pm(i)/rhop(i)
         TEp(i)=2.*TEp(i)-TEo(i)
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i)) 
      end do


      if(iDBC.eq.1)then	!Hughes & Grahams 2009 correction 
           do i=1,nb
            rhop(i) = 2.*rhop(i) - rhoo(i)
            if (rhop(i).lt.rho0) rhop(i)=rho0
            pVol(i) = pm(i)/rhop(i)
            TEp(i)=2.*TEp(i)-TEo(i)
            call equation_of_state(rhop(i),TEp(i),p(i),cs(i)) 
           enddo
      endif 

     
      do i=nbp1,np
         xp(i)   = 2.*xp(i) - xo(i)
         zp(i)   = 2.*zp(i) - zo(i)

         up(i)   = 2.*up(i) - uo(i)
         wp(i)   = 2.*wp(i) - wo(i)
         
      enddo

       
      call movingObjects(1)

      !--- Diagnostic Print out ---
      if(0.eq.i_rank_check.and.itime.ge.itime_check)then
          i = i_ParticleCheck  
          print*
          print*,'Accn Check, AFTER CORRECT, CORRECTOR'
          print*,'rank_id, itime ',rank_id, itime
          print*,'i',i
          print*,'xp(i)',xp(i)
          print*,'zp(i)',zp(i)
          print*,'up(i)',up(i)
          print*,'wp(i)',wp(i)
          print*,'p(i)',p(i)
          print*,'rhop(i)',rhop(i)
          print*,'udot(i)',udot(i)
          print*,'wdot(i)',wdot(i)
          print*,'rdot(i)',rdot(i)
          print*,'dudx_CSPH(i)',dudx_CSPH(i)
          print*,'dudz_CSPH(i)',dudz_CSPH(i)
          print*,'dwdx_CSPH(i)',dwdx_CSPH(i)
          print*,'dwdz_CSPH(i)',dwdz_CSPH(i)
          print*,'drhodx_CSPH(i)',drhodx_CSPH(i)
          print*,'drhodz_CSPH(i)',drhodz_CSPH(i)
          !print*,'sum_wab(i) ',sum_wab(i)
          !print*,'sum_dWdx(i)',sum_dWdx(i)
          !print*,'sum_dWdy(i)',sum_dWdy(i)
          !print*,'rhop_sum(i)',rhop_sum(i)
          print*,'ddt_p, ddt_c, dt',ddt_p, ddt_c, dt
          if(itime.ge.itime_check+0)stop
      end if

      end

