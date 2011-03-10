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


      subroutine rigid_body_motion(i_midPoint)

      include 'common.2D'
      
      double precision dxp_dble1, dxp_dble2 

      real(kind=kind_XZ) Box_XC_SOI, Box_ZC_SOI
           
      itime_check = 10e8     !0    !10e8     !0    !10e8     !
      num_FB_check = 10
      i_check = 1362
      !- Loop over bodies individually -
      do i_num_FB = 1,num_FB 
       
c        if(i_midPoint.ge.0)then   !if(i_midPoint.eq.0)then
          bigU_old(i_num_FB)      = BigU(i_num_FB)  !
          bigW_old(i_num_FB)      = BigW(i_num_FB)

          !bigOmegaY_old(i_num_FB) = bigOmegaY(i_num_FB)
          Box_XC_old(i_num_FB)    = Box_XC(i_num_FB)
          Box_ZC_old(i_num_FB)    = Box_ZC(i_num_FB)

c        endif
        bigOmegaY_old(i_num_FB) = bigOmegaY(i_num_FB)
        


        !Diagnostic info
        !i = nb_FB(i_num_FB) - 2
        !xo_temp = xp(i)
        !zo_temp = zp(i)
        !uo_temp = up(i)
        !wo_temp = wp(i)


        
        if(i_num_FB.eq.1)then
          i_start = nbfm + 1
        else
          i_start = nb_FB(i_num_FB-1) + 1
        end if


        OmegaYdot(i_num_FB) = 0.0

                
        bigUdot_init      = X_Friction(i_num_FB)  !Friction forces summed so far (not per unit mass)
        bigWdot_init      = Z_Friction(i_num_FB)  !Friction forces summed so far (not per unit mass)

        
        !Uncomment these if the body forces are NOT calculated in monBC
        !bigUdot(i_num_FB) = X_Friction(i_num_FB)  !Friction forces summed so far (not per unit mass)
        !bigWdot(i_num_FB) = Z_Friction(i_num_FB)  !Friction forces summed so far (not per unit mass)
        !
        !BodyForce_x_total = X_Friction(i_num_FB)
        !BigU_old        = BigU(i_num_FB)*i_midPoint
        !BigW_old        = BigW(i_num_FB)*i_midPoint
        !
        !bigOmegaY_old   = bigOmegaY(i_num_FB)*i_midPoint
        !Box_XC_old      = Box_XC(i_num_FB)*i_midPoint
        !Box_ZC_old      = Box_ZC(i_num_FB)*i_midPoint


        do i = i_start,nb_FB(i_num_FB)
        
          bigUdot(i_num_FB) = bigUdot(i_num_FB) + ax(i)   !This accn is not per unit mass
          bigWdot(i_num_FB) = bigWdot(i_num_FB) + az(i)   !This accn is not per unit mass



          OmegaYdot(i_num_FB) = OmegaYdot(i_num_FB)
     &          + (zp_minus_R0(i)*ax(i) - xp_minus_R0(i)*az(i))
     

     
        end do

        !-- Wave force on one side of block --
        BodyForce_x = 0.0
        BodyForce_z = 0.0

        do i = i_start,nb_FB(i_num_FB)        
          if(xp(i).lt.Box_XC(i_num_FB))then
            BodyForce_x = BodyForce_x + ax(i)   !This accn is not per unit mass
          endif
          !BodyForce_x = BodyForce_x + ax(i)   !This accn is not per unit mass
          BodyForce_z = BodyForce_z + az(i)   !This accn is not per unit mass
          BodyForce_x_total = BodyForce_x_total + ax(i)   !This accn is not per unit mass
        end do
        BodyForce_x = BodyForce_x + bigMass(i_num_FB)*grx
        BodyForce_z = BodyForce_z + bigMass(i_num_FB)*grz
 
        BodyForce_x_temp = BodyForce_x
        BodyForce_z_temp = BodyForce_z

        BodyForce_x_old1 = BodyForce_x_temp
        BodyForce_z_old1 = BodyForce_z_temp
        

        
        !- Include Body Forces -
        bigUdot(i_num_FB) = bigUdot(i_num_FB) + bigMass(i_num_FB)*grx
        bigWdot(i_num_FB) = bigWdot(i_num_FB) + bigMass(i_num_FB)*grz

        
        !- Evaluate non-friction forces -
        if(nb_inFriction(i_num_FB).gt.0)then
          X_nonFriction(i_num_FB) = 
     &         (bigUdot(i_num_FB) - X_Friction(i_num_FB))                     !Total non-friction forces
     &         /(nb_inFriction(i_num_FB))                   !per particle involved in Friction
c         Y_nonFriction(i_num_FB) = 
c     &        (bigVdot(i_num_FB) - Y_Friction(i_num_FB))                     !Total non-friction forces
c     &        /(nb_inFriction(i_num_FB))                   !per particle involved in Friction
          Z_nonFriction(i_num_FB) = 
     &         (bigWdot(i_num_FB) - Z_Friction(i_num_FB))                     !Total non-friction forces
     &         /(nb_inFriction(i_num_FB))                   !per particle involved in Friction
        else
          X_nonFriction(i_num_FB) = 
     &         (bigUdot(i_num_FB) + bigMass(i_num_FB)*grx 
     &          - X_Friction(i_num_FB))                     !Total non-friction forces
c         Y_nonFriction(i_num_FB) = 
c     &        (bigVdot(i_num_FB) + bigMass(i_num_FB)*gry 
c     &         - Y_Friction(i_num_FB))                     !Total non-friction forces
          Z_nonFriction(i_num_FB) = 
     &         (bigWdot(i_num_FB) + bigMass(i_num_FB)*grz 
     &          - Z_Friction(i_num_FB))                     !Total non-friction forces
        endif
c        if(i_num_FB.eq.-num_FB_check)then
c          print*
c          print*,'rigid_body_motion '
c          print*,'itime, i_num_FB ',itime, i_num_FB 
c          print*,'Box_XC, Box_ZC ',Box_XC(i_num_FB),Box_ZC(i_num_FB)
c          print*,'bigUdot(i_num_FB) ',bigUdot(i_num_FB)
c          print*,'bigWdot(i_num_FB) ',bigWdot(i_num_FB)
c          print*,'bigMass(i_num_FB) ',bigMass(i_num_FB)
c          print*,'grx, grz ',grx, grz
c          print*,'X_Friction(i_num_FB)    ',X_Friction(i_num_FB)
c          print*,'Z_Friction(i_num_FB)    ',Z_Friction(i_num_FB)
c          print*,'X_nonFriction(i_num_FB) ',X_nonFriction(i_num_FB)
c          print*,'Z_nonFriction(i_num_FB) ',Z_nonFriction(i_num_FB)
c          print*,'nb_inFriction(i_num_FB) ',nb_inFriction(i_num_FB)
c          read(*,*)
c          print*
c        endif
        
        !- Forces per unit mass -
        bigUdot(i_num_FB) = bigUdot(i_num_FB)/bigMass(i_num_FB) !+ grx
        bigWdot(i_num_FB) = bigWdot(i_num_FB)/bigMass(i_num_FB) !+ grz

        OmegaYdot(i_num_FB) = OmegaYdot(i_num_FB)/bigInertiaYY(i_num_FB)
        
        
        Box_XC(i_num_FB) = Box_XC_old(i_num_FB) + dt2*bigU(i_num_FB)
        Box_ZC(i_num_FB) = Box_ZC_old(i_num_FB) + dt2*bigW(i_num_FB)

        bigU(i_num_FB) = bigU_old(i_num_FB) + dt2*bigUdot(i_num_FB)
        bigW(i_num_FB) = bigW_old(i_num_FB) + dt2*bigWdot(i_num_FB)



        bigOmegaY(i_num_FB) = bigOmegaY_old(i_num_FB)
     &                    + dt2*OmegaYdot(i_num_FB)
        


        !- Account for periodicity
        if(i_periodicOBs(1).eq.1)then  !Replace Box Centre on other side of domain
          if(Box_XC(i_num_FB).gt.xmax_container)then
            dxp_dble1 = dble(Box_XC(i_num_FB)) - xmax_container_double 
            Box_XC(i_num_FB) = real(dxp_dble1 +xmin_container_double)
          else if(Box_XC(i_num_FB).lt.xmin_container)then
            dxp_dble1 = xmin_container_double - dble(Box_XC(i_num_FB))   
            Box_XC(i_num_FB) = real(xmax_container_double - dxp_dble1)
          end if
        endif

        do i = i_start,nb_FB(i_num_FB) 

         xdot(i) = up(i)
         zdot(i) = wp(i)


         xp(i) = xp(i) + dt2*xdot(i)
         zp(i) = zp(i) + dt2*zdot(i)


         xp_minus_R0(i) = xp(i) - Box_XC(i_num_FB)
         zp_minus_R0(i) = zp(i) - Box_ZC(i_num_FB)


         !- Account for periodicity
         !THIS NEEDS TO BE DONE MORE ROBUSTLY!!
         if(i_periodicOBs(1).eq.1)then
           if(xp_minus_R0(i).gt.(0.5*vlx))then
             dxp_dble1 = xmax_container_double-dble(xp(i))
             dxp_dble2 = dble(Box_XC(i_num_FB))-xmin_container_double
             xp_minus_R0(i) = -real(dxp_dble1+dxp_dble2)
           end if
           if(-xp_minus_R0(i).gt.(0.5*vlx))then
             dxp_dble1 = dble(xp(i))-xmin_container_double
             dxp_dble2 = xmax_container_double-dble(Box_XC(i_num_FB))
             xp_minus_R0(i) =  real(dxp_dble1+dxp_dble2)
           end if
         endif

         up(i) = bigU(i_num_FB) + bigOmegaY(i_num_FB)*zp_minus_R0(i)
         wp(i) = bigW(i_num_FB) - bigOmegaY(i_num_FB)*xp_minus_R0(i)


         

        enddo
        
      end do   !End of:  do i_num_FB = 1,num_FB

      if(iBC.eq.1)then
        call updateNormals_2D(nbfm + 1,nb)        
      endif

      end subroutine
