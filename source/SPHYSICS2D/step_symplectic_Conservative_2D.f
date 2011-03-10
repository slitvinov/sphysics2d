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

c     Symplectic Method

c  ...  compute acceleration
c
c      !-- Debug Code------
c      itime_check = 10e8   !0    !10e8   !0    !
c      i_ParticleCheck = 1722   !1101   !1284    !nb + 1      

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
      
c      !-- Debug Code------
c      if(itime.ge.itime_check)then
c        print*
c        print*,'AFTER CORRECT, PREDICTOR, before movingOBJECTS'
c        print*,'rank_id, itime, dt ',rank_id, itime, dt
c        print*,'ddt_p, dt, dt2 ',ddt_p, dt, dt2
c        i = i_ParticleCheck  
c        print*,'i = ',i
c        print*,'xp, zp  ',xp(i),zp(i)
c        print*,'up, wp  ',up(i),wp(i)
c        print*,'rhop, p,  pVol, pm ',rhop(i),p(i) ,pVol(i),pm(i)
c        print*,'rhoo, po, Volo, pm ',rhoo(i),po(i),Volo(i),pm(i)
c        print*,'VolTotalE   , VolTotalEo ',
c     &          VolTotalE(i), VolTotalEo(i) 
c        print*,'aTE(i) ',aTE(i)
c        print*,'TEp, TotalE ',TEp(i),TotalE(i)
c        print*,'ar(i), aVol(i) ',ar(i), aVol(i)
c        print*,'Volrhoo(i), Volrho(i) ',Volrhoo(i), Volrho(i) 
c        print*,'ax(i), Volrhou(i) ',ax(i), Volrhou(i)
cc       print*,'ay(i), Volrhov(i) '       
c        print*,'az(i), Volrhow(i) ',az(i), Volrhow(i)
c        print*,'pVol(i), Volo(i) ',pVol(i), Volo(i)
c        print*,'dudx_CSPH(i)',dudx_CSPH(i)
c        print*,'dudz_CSPH(i)',dudz_CSPH(i)
c        print*,'dwdx_CSPH(i)',dwdx_CSPH(i)
c        print*,'dwdz_CSPH(i)',dwdz_CSPH(i)
c        print*,'drhodx_CSPH(i)',drhodx_CSPH(i)
cc       print*,'drhody_CSPH(i)',drhody_CSPH(i)
c        print*,'drhodz_CSPH(i)',drhodz_CSPH(i)
c        print*,'sum_wab(i) ',sum_wab(i)
c        write(*,859)'ax(i) ',ax(i)
cc       write(*,859)'ay(i) '
c        write(*,859)'az(i) ',az(i)
c        if(i.gt.nbfm.and.i.le.nb)then
c          i_num_FB = i_FB_Pointer_Info(i)
c          print*,'i, i_num_FB ',i, i_num_FB
c          print*,'Box_XC, Box_ZC         ',
c     &            Box_XC(i_num_FB), Box_ZC(i_num_FB)
c          print*,'bigU, bigW, bigOmegaY ',
c     &            bigU(i_num_FB), bigW(i_num_FB), bigOmegaY(i_num_FB) 
c          print*,'bigUdot, bigWdot ',
c     &            bigUdot(i_num_FB), bigWdot(i_num_FB)
c          print*,'OmegaYdot ',
c     &            OmegaYdot(i_num_FB) 
c          print*,'X_Friction, Z_Friction ',
c     &            X_Friction(i_num_FB), Z_Friction(i_num_FB)
c          print*,'X_nonFriction, Z_nonFriction ',
c     &            X_nonFriction(i_num_FB), Z_nonFriction(i_num_FB)
c          print*,'nb_inFriction ',
c     &            nb_inFriction(i_num_FB) 
c          print*,'Box_XC_old, Box_ZC_old   ',
c     &            Box_XC_old(i_num_FB), Box_ZC_old(i_num_FB)
c          print*,'bigU_old, bigW_old       ',
c     &            bigU_old(i_num_FB), bigW_old(i_num_FB)
c          print*,'xp(i)  , zp(i)  ',xp(i), zp(i)   
c          print*,'up(i)  , wp(i)  ',up(i), wp(i)  
c          print*,'ax(i)  , udot(i)',ax(i), udot(i)
c          print*,'az(i)  , wdot(i)',az(i), wdot(i)
c          print*,'xp_minus_R0(i), zp_minus_R0(i) ',
c     &            xp_minus_R0(i), zp_minus_R0(i)
c          print*,'- OmegaY*zp_minus_R0(i) ',
c     &            - bigOmegaY(i_num_FB)*zp_minus_R0(i)
c          print*,'  OmegaY*xp_minus_R0(i) ',
c     &              bigOmegaY(i_num_FB)*xp_minus_R0(i)
c        endif
c      endif

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
      do i=nstart,np
         epsilon_Voldot = aVol(i)/pVol(i)
         pVol(i) = Volo(i)*
     &            (2.0 + epsilon_Voldot*dt)/(2.0 - epsilon_Voldot*dt)
         Volrho(i) = Volrhoo(i)*
     &     (2.0*Volrho(i) + ar(i)*dt)/(2.0*Volrho(i) - ar(i)*dt)   !This one ??!!
         rhop(i) = Volrho(i) / pVol(i)
      end do
      
      do i=nbp1,np
         Volrhou(i)   = Volrhouo(i) + dt*ax(i)
         Volrhow(i)   = Volrhowo(i) + dt*az(i)
         up(i) = Volrhou(i) / Volrho(i)
         wp(i) = Volrhow(i) / Volrho(i)
      
         xp(i)   = xo(i) + dt2*(uo(i) + up(i))
         zp(i)   = zo(i) + dt2*(wo(i) + wp(i))
         
      enddo

      do i=nstart,np
c  ...   Compute thermal energy
         VolTotalE(i) = VolTotalE(i) + dt2*aTE(i)
         TotalE(i) = VolTotalE(i) / pVol(i)
         TEp(i) = TotalE(i)/rhop(i) - 0.5*(up(i)*up(i) + wp(i)*wp(i))
c  ...   compute pressure using the equation of state
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
      enddo

      call movingObjects(1)

      if (ivar_dt.eq.1) then
          call variable_time_step 
         ddt_c=dt_new
         dt=min(ddt_p,ddt_c)  !Time step to be used in next loop
      endif
       
c      if(itime.ge.itime_check)then
c        print*
c        print*,'Accn Check, AFTER CORRECT, CORRECTOR'
c        print*,'itime ',itime
c        print*,'ddt_p, ddt_c, dt, dt2 ',ddt_p, ddt_c, dt, dt2
c        i = i_ParticleCheck  
c        print*,'i ',i
c        print*,'xp, zp  ',xp(i),zp(i)
c        print*,'up, wp  ',up(i),wp(i)
c        print*,'rhop, p,  pVol, pm ',rhop(i),p(i) ,pVol(i),pm(i)
c        print*,'rhoo, po, Volo, pm ',rhoo(i),po(i),Volo(i),pm(i)
c        print*,'VolTotalE   , VolTotalEo ',
c     &          VolTotalE(i), VolTotalEo(i) 
c        print*,'aTE(i) ',aTE(i)
c        print*,'TEp, TotalE ',TEp(i),TotalE(i)
c        print*,'ar(i), aVol(i) ',ar(i), aVol(i)
c        print*,'Volrhoo(i), Volrho(i) ',Volrhoo(i), Volrho(i) 
c        print*,'ax(i), Volrhou(i) ',ax(i), Volrhou(i)
c        print*,'az(i), Volrhow(i) ',az(i), Volrhow(i)
c        print*,'pVol(i), Volo(i) ',pVol(i), Volo(i)
c        print*,'dudx_CSPH(i)',dudx_CSPH(i)
c        print*,'dudz_CSPH(i)',dudz_CSPH(i)
c        print*,'dwdx_CSPH(i)',dwdx_CSPH(i)
c        print*,'dwdz_CSPH(i)',dwdz_CSPH(i)
c        print*,'drhodx_CSPH(i)',drhodx_CSPH(i)
c        print*,'drhodz_CSPH(i)',drhodz_CSPH(i)
c        print*,'sum_wab(i) ',sum_wab(i)
c        write(*,859)'ax(i) ',ax(i)
c        write(*,859)'az(i) ',az(i)
c        if(i.gt.nbfm.and.i.le.nb)then
c          i_num_FB = i_FB_Pointer_Info(i)
c          print*,'i, i_num_FB ',i, i_num_FB
c          print*,'Box_XC, Box_ZC         ',
c     &            Box_XC(i_num_FB), Box_ZC(i_num_FB)
c          print*,'bigU, bigW, bigOmegaY ',
c     &            bigU(i_num_FB), bigW(i_num_FB), bigOmegaY(i_num_FB) 
c          print*,'bigUdot, bigWdot ',
c     &            bigUdot(i_num_FB), bigWdot(i_num_FB)
c          print*,'OmegaYdot ',
c     &            OmegaYdot(i_num_FB) 
c          print*,'X_Friction, Z_Friction ',
c     &            X_Friction(i_num_FB), Z_Friction(i_num_FB)
c          print*,'X_nonFriction, Z_nonFriction ',
c     &            X_nonFriction(i_num_FB), Z_nonFriction(i_num_FB)
c          print*,'nb_inFriction ',
c     &            nb_inFriction(i_num_FB) 
c          print*,'Box_XC_old, Box_ZC_old   ',
c     &            Box_XC_old(i_num_FB), Box_ZC_old(i_num_FB)
c          print*,'bigU_old, bigW_old       ',
c     &            bigU_old(i_num_FB), bigW_old(i_num_FB)
c          print*,'xp(i)  , zp(i)  ',xp(i), zp(i)   
c          print*,'up(i)  , wp(i)  ',up(i), wp(i)  
c          print*,'ax(i)  , udot(i)',ax(i), udot(i)
c          print*,'az(i)  , wdot(i)',az(i), wdot(i)
c          print*,'xp_minus_R0(i), zp_minus_R0(i) ',
c     &            xp_minus_R0(i), zp_minus_R0(i)
c          print*,'- OmegaY*zp_minus_R0(i) ',
c     &            - bigOmegaY(i_num_FB)*zp_minus_R0(i)
c          print*,'  OmegaY*xp_minus_R0(i) ',
c     &              bigOmegaY(i_num_FB)*xp_minus_R0(i)
c        endif
c        if(itime.eq.0)stop
cc        stop
c      end if
c859     format(a11,' ',e24.16)        

      end
c
