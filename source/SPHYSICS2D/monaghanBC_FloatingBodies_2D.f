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


      subroutine monaghanBC(iip,jjp,drx,drz,dot,
     &                      duxp,duzp,
     &                      fxbp,fzbp)

      include 'common.2D'
      id_num_FB_i = 0
      id_num_FB_j = 0
c     Boundary condition JWPCOE, 32
c
c     iip is a WaterParticle (WP)
c     jjp is a BoundaryParticle (BP)
c
c     drxp = xp(WP) - xp(BedP)
c     drzp = zp(WP) - zp(BedP)


c      Projection of WP position onto normal
       pdot3 = drx*xnb(jjp)+drz*znb(jjp)

c       if(pdot3.gt.0.0)then  !ensures that WP is on correct side of boundary
c
c         xpp = magnitude(vec(dr) x vec(n))   cross-product
c
c         xpp=(drxp*ynb(jjp)-dryp*xnb(jjp))


c        Projection of WP on plane of boundary
c
c        r(BP->WP) = [drxp,dryp,drzp]
c
c        r(BP->WP projection) = r(BP->WP) -[r(BP->WP).n/(|n|.|n|)]n
c
c        normal to plane = n
c        magnitude of normal = |n| = 1

         drx_proj = drx - (pdot3)*xnb(jjp)
         drz_proj = drz - (pdot3)*znb(jjp)


c        Projection of plane-projected vector r(BP->WP projection) onto tangent

         xpp = drx_proj*xtb(jjp)+drz_proj*ztb(jjp)


         if(xpp.gt.0.0)then
           deltappt = deltaptb(jjp,2)
         else
           deltappt = deltaptb(jjp,1)
         end if
c
c
c
c
c      
c
         Prx = abs(xpp)/deltappt    
         
         facp = 0.0
         f_tb = 0.0
         f_fb = 0.0
         chi_X = 0.0

         Q_fn = 0.0
         i_CornerParticle = 0

         !normalized distance from wall
         if(pdot3.gt.0.0)then ! particle inside
             q=pdot3 * one_over_2h           !normalized distance from wall
         else if(pdot3.lt.0.0)then !Particle ouside 
             q=0.2*abs(pdot3) * one_over_2h  !normalized distance from wall
         else !Particle on the boundary
             q=h*(1.0e-3) * one_over_2h      !correcting for zero distance from wall
         end if

c         if(q.lt.1.0.and.pdot3.gt.0.0.and.
c     &                abs(xpp).lt.deltappt)then
         if(q.lt.1.0.and.abs(xpp).lt.deltappt
     &              .and.((iip.gt.nb.and.pdot3.gt.0.0).or.
     &                    (iip.lt.nbp1)))then
     
             u_normal = duxp*xnb(jjp) + duzp*znb(jjp)
c            u_normal < 0 if particles are approaching
c            u_normal > 0 if particles are moving apart
      
             !-- Corner Particle Identification --

             i_minus1 = iBP_Pointer_Info(jjp,3)     !Absolute index of i-1 neighbour BP
             i_plus1  = iBP_Pointer_Info(jjp,4)     !Absolute index of i+1 neighbour BP
c
c
             
             !c0 = coeff.Vmax (Vmax = sqrt(g.h_SWL); coef = beta_coef.10, beta_coef2 = beta_coef**2
             !beta_coef2 = Factor that B is changed by 
             if(iip.gt.nb)then  !Normal FP-BP Interaction
               if(u_normal.gt.0.0)then
                 App=(0.01*cs(iip)*cs(iip)*one_over_beta_coef2
     &                 -cs(iip)*u_normal*one_over_beta_coef)*one_over_h
               else
                 App=0.01*(cs(iip)*cs(iip)*one_over_beta_coef2)
     &                   *one_over_h    
               endif
             else !BP-BP interactions
               cs_i = cs0
               if(u_normal.gt.0.0)then    !Particles receding:  u_normal > 0
                 cs2_diff = 0.5*abs(cs_i*cs_i)   !
                 App=(0.01*cs2_diff*one_over_beta_coef2
     &                 -cs_i*u_normal*one_over_beta_coef)*one_over_h

               else                       !Particles approaching:  u_normal < 0
                 !cs2_diff = 0.5*abs(cs_i*cs_i)   !Set this for elastic bouncing
                 cs2_diff = 0.0   !Set this for elastic bouncing
                 repulsion_factor = 10.0
                 !if(iip.gt.nbfm.and.jjp.gt.nbfm)then
                 !  repulsion_factor = 2.0
                 !endif
                 App=(0.01*cs2_diff*one_over_beta_coef2
     &                 -repulsion_factor*cs_i*u_normal
     &                 *one_over_beta_coef)*one_over_h

               endif
               !App = App*0.5
               if(q.lt.0.2.and.
     &                   (jjp.ne.i_minus1.and.
     &                    jjp.ne.i_plus1))then


                 if(q.lt.0.2)then !Boundary Particle has crossed boundary 
                   q = h*(1.0e-3) * one_over_2h 
                 endif
                 q = q*q
               end if
             endif

	
c          --- Depth Correction ---
             eps_min = 0.02
             epsilon_Z = abs(1.0 - (zp(jjp)/h_SWL))
             dyn_head =(up(jjp)*up(jjp) + wp(jjp)*wp(jjp))
     &                  /(2.0*abs(grz)*h_SWL)
c     &                  /(2.0*9.81*h_SWL)   !1-D
           if(iip.gt.nb)then  !WP-BP
             if(epsilon_Z.gt.1.0-eps_min)then !Close to the bottom (No correction )
                epsilon_Z = 1.0
             elseif(zp(iip).gt.h_SWL)then !over still water level          
                epsilon_Z = eps_min + dyn_head
             else                 !under still water level & far from bottom
                epsilon_Z = eps_min + epsilon_Z + dyn_head
             endif
           else   !BP-BP
             epsilon_Z = 1.0 + dyn_head
           endif

c          --- Velocity Correction ---
             if(u_normal.gt.0.0)then
                epsilon_dyn = 0.0
             else 
                epsilon_dyn = 40.*abs(u_normal)/cs0     !factor 20 needed 
			                                  ! to keep scheme stable 											      ! with Monaghan BC wavemaker
             endif

             if(epsilon_dyn.gt.1.0)then
                epsilon_dyn = 1.0
             endif

             Appc = App*(epsilon_Z + epsilon_dyn) !Force correction
             Q_fn = (1.0-q)/sqrt(q)
             
             facp=0.5*(1.+cos(pi*abs(xpp)/deltappt))*
     1            (Appc*Q_fn)   !(Appc*(1.0-q)/sqrt(q))
     2            *vnorm_mass


             fxbp = xnb(jjp)*facp
             fzbp = znb(jjp)*facp   


             if(iip.lt.nbp1)then
               fmax = Appc*Q_fn*vnorm_mass
               chi_X = (1.0 - Prx)

               q_r = q
               xtb_r = xtb(jjp)
               ztb_r = ztb(jjp)

               abs_t_dot_t = abs(xtb(iip)*xtb_r 
     &                         + ztb(iip)*ztb_r)
c
c
               !-- Friction force between rigid objects --

               id_num_FB_i = i_FB_Pointer_Info(iip)
               pmu = friction_coeff(id_num_FB_i)
               if(abs_t_dot_t.lt.0.6)then  !Activates forces only between sliding faces
                 pmu = 0.0
                 fxbp = 0.0
                 fzbp = 0.0

                 facp_F = 0.0
               else
                 facp_F = fmax
               endif
               !if(jjp.gt.nbfm)pmu = 0.0  !No friction between moving blocks   !<--------------------
               !u_parallel = duxp*xtb(jjp) + duzp*ztb(jjp)
               u_parallel = duxp*xtb_r + duzp*ztb_r

               if(u_parallel.gt.0.0)then
                 sign_u_parallel = 1.0
               elseif(u_parallel.lt.0.0)then
                 sign_u_parallel = -1.0
               else
                 sign_u_parallel = 0.0
               endif
               !f_fb = - sign_u_parallel*pmu*facp
               
               F_Static_i = (X_nonFriction(id_num_FB_i)*xtb_r
     &                     + Z_nonFriction(id_num_FB_i)*ztb_r)


               F_Static_i_abs = abs(F_Static_i)

               if(F_Static_i.gt.0.0)then
                 sign_F_Static_i = 1.0
               elseif(F_Static_i.lt.0.0)then
                 sign_F_Static_i = -1.0
               else
                 sign_F_Static_i = 0.0
               endif
               if(sign_u_parallel.eq.sign_F_Static_i)then
                 !F_Static_abs = F_Static_abs * 1.0005   !An extra 2.5% for cases where the block inches forward
               endif
               
               if(iip.gt.nbfm)then  !iip is a floating object particle
                 id_num_FB_i = i_FB_Pointer_Info(iip)
                 F_Dynamic = pmu*facp_F*pm(iip)
                 u_parallel_ratio = abs(u_parallel)
     &                             *one_over_u_parallel_max
                 if(u_parallel_ratio.lt.1.0)then
                   !- Static Friction Force -
                   f_temp = - sign_u_parallel*chi_X
     &             *min(F_Dynamic,F_Static_i_abs)  !ensures max friction force is <= F_Dynamic
                   F_fb_x_i =   f_temp*xtb_r
                   F_fb_z_i =   f_temp*ztb_r

                 elseif(u_parallel_ratio.lt.2.0)then
                   !- Static-to-Dynamic Friction Force F_SD -
                   F_SD = F_Static_i_abs + (F_Dynamic - F_Static_i_abs)
     &                           *(u_parallel_ratio - 1.0)  !linear interpolation
                   f_temp = - sign_u_parallel*chi_X
     &             *min(F_Dynamic,F_SD)  !ensures max friction force is <= F_Dynamic
                   F_fb_x_i =   f_temp*xtb_r
                   F_fb_z_i =   f_temp*ztb_r

                 else
                   !- Dynamic Friction Force -
                   f_temp = - sign_u_parallel*chi_X*F_Dynamic
                   F_fb_x_i =   f_temp*xtb_r
                   F_fb_z_i =   f_temp*ztb_r

                 endif
                 X_Friction(id_num_FB_i) = X_Friction(id_num_FB_i)
     &                                        + F_fb_x_i
c                Y_Friction(id_num_FB_i) = Y_Friction(id_num_FB_i)
c    &                                        + F_fb_y_i
                 Z_Friction(id_num_FB_i) = Z_Friction(id_num_FB_i)
     &                                        + F_fb_z_i
                 fxbp = fxbp + F_fb_x_i/pm(iip)  !per unit mass force
                 fzbp = fzbp + F_fb_z_i/pm(iip)  !per unit mass force

                 if(chi_X.gt.0.5.and.abs_t_dot_t.gt.0.6)then
                   nb_inFriction(id_num_FB_i) = 
     &             nb_inFriction(id_num_FB_i) + 1
                 endif
c                      if(iip.eq.-1135)then
c                       print*,'I,j,xp,zp, num_FB, nb_inFriction#',
c     &                 iip,jjp,xp(iip),zp(iip),xp(jjp),zp(jjp),
c     &                 id_num_FB_i,nb_inFriction(id_num_FB_i),
c     &                 duxp,duzp,u_parallel,
c     &                 u_parallel_ratio,xtb_r,ztb_r,
c     &                 sign_u_parallel,chi_X,F_Dynamic,F_Static_i_abs,
c     &                 pmu,facp_F,pm(iip),abs_t_dot_t,
c     &                 f_temp,F_fb_x_i,F_fb_z_i,
c     &                 X_Friction(id_num_FB_i),Z_Friction(id_num_FB_i)
c                     endif
                !- This section needs to be updated for multiple floating bodies sliding past each other -
                 if(jjp.gt.nbfm)then    !jjp is also a floating object particle
                   id_num_FB_j = i_FB_Pointer_Info(jjp)
                   F_Static_j = (X_nonFriction(id_num_FB_j)*xtb_r
     &                         + Z_nonFriction(id_num_FB_j)*ztb_r)

                   F_Static_j_abs = abs(F_Static_j)
                   F_fb_x_j = F_fb_x_i
                   F_fb_z_j = F_fb_z_i

                   X_Friction(id_num_FB_j) = X_Friction(id_num_FB_j)
     &                                          - F_fb_x_j
c                  Y_Friction(id_num_FB_j) = Y_Friction(id_num_FB_j)
c    &                                          - F_fb_y_j
                   Z_Friction(id_num_FB_j) = Z_Friction(id_num_FB_j)
     &                                          - F_fb_z_j
                   if(chi_X.gt.0.5.and.abs_t_dot_t.gt.0.6)then
                     nb_inFriction(id_num_FB_j) = 
     &               nb_inFriction(id_num_FB_j) + 1
c                     if(jjp.eq.1135)then
c                       print*,'i,J,xp,zp, num_FB, nb_inFriction#',
c     &                 iip,jjp,xp(iip),zp(iip),xp(jjp),zp(jjp),
c     &                 id_num_FB_j,nb_inFriction(id_num_FB_j)
c                     endif
                   endif
                 endif   !End of:    if(jjp.gt.nbfm.and.iip.lt.nbp1)then
               endif   !if(iip.gt.nbfm)then
             endif   !End of:    if(1.eq.2.and.iip.lt.nbp1)then
         else
             fxbp = 0.0
             fzbp = 0.0

         end if

c      if(1.eq.2.and.
c     &  itime.ge.0.and.   !itime.ge.11100.and.    !itime.ge.0.and.   ! 
c     &  jjp.le.nb.and.(iip.eq.1093.or.iip.eq.1103))then
cc        if(q.lt.1.0.and.   !pdot3.lt.0.1*h.and.
cc     &      abs(xpp).lt.deltappt)then
c         if(jjp.gt.nbfm)then
c           j_num_FB = i_FB_Pointer_Info(jjp) 
c         else
c           j_num_FB = 0 
c         endif
c         print*
c         print*,'rank_id, itime ',rank_id, itime
c         print*,'iip  ',iip
c         print*
c         print*,'i_num_FB ',i_num_FB
c         print*,'xp,zp',xp(iip),zp(iip)
c         print*,'up,wp',up(iip),wp(iip)
c         print*,'jjp ',jjp
c         print*
c         print*,'j_num_FB ',j_num_FB
c         print*,'xp,zp',xp(jjp),zp(jjp)
c         print*,'up,wp',up(jjp),wp(jjp)
c         print*,'xnb,znb',xnb(jjp),znb(jjp)
c         print*,'xtb,ztb',xtb(jjp),ztb(jjp)
c         print*,'drx ,drz ',drx ,drz
c         print*,'itime ',itime
c         print*,'pdot3',pdot3
c         print*,'one_over_2h ',one_over_2h
c         print*,'q',q
c         print*,'xpp',xpp
c         print*,'deltappt',deltappt
c         print*,'rhop(iip) ',rhop(iip)
c         print*,'cs(iip)',cs(iip)
c         print*,'epsilon_Z',epsilon_Z
c         print*,'epsilon_dyn',epsilon_dyn
c         print*,'duxp, duzp, u_normal ', duxp, duzp, u_normal
c         print*,'one_over_beta_coef2, one_over_beta_coef, one_over_h ',
c     &           one_over_beta_coef2, one_over_beta_coef, one_over_h
c         print*,'cs0, cs_i, cs2_diff ',cs0, cs_i, cs2_diff
c         print*,'App ',App
c         print*,'Appc',Appc
c         print*,'facp',facp
c         print*,'App*(1-q)/sqrt(q)',App*(1.0-q)/sqrt(q)
c         print*,'facp ',facp
c         print*,'fxbp ',fxbp
c         print*,'fzbp ',fzbp
c         print*
c         !read(*,*)
cc        end if
c        !stop
c      end if
c      if(1.eq.2.and.
c     &  itime.ge.2.and.   !itime.ge.11100.and.    !itime.ge.0.and.   !
c     &  jjp.le.nb.and.(iip.eq.1362.or.jjp.eq.1362))then
c        if(q.lt.1.0.and.   !pdot3.lt.0.1*h.and.
c     &      abs(xpp).lt.deltappt)then
c          if(id_num_FB_i.eq.10.or.id_num_FB_j.eq.10)then
c            print*
c            print*,'rank_id, itime ',0, itime
c            print*,'iip                ',iip
c            print*,'i_ParticleType_iip '
c            print*,'id_num_FB_i ',id_num_FB_i
c            print*,'xp,zp',xp(iip),zp(iip)
c            print*,'up,wp',up(iip),wp(iip)
c            print*,'jjp               ',jjp
c            print*,'i_ParticleType_jjp '
c            print*,'id_num_FB_j ',id_num_FB_j
c            print*,'xp,zp',xp(jjp),zp(jjp)
c            print*,'up,wp',up(jjp),wp(jjp)
c            print*,'xnb,znb',xnb(jjp),znb(jjp)
c            print*,'xtb,ztb',xtb(jjp),ztb(jjp)
c            print*,'drx ,drz ',drx ,drz    
c            print*,'xpp, deltappt',xpp, deltappt
c            print*,'i_plus1, i_minus1 ',i_plus1, i_minus1
c            print*,'pdot3, q ',pdot3, q
c            print*,'one_over_2h ',one_over_2h
c            print*,'Appc, facp ',Appc, facp
c            print*,'chi_X, abs_t_dot_t',chi_X, abs_t_dot_t
c            print*,'u_parallel, u_parallel_ratio, sign_u_parallel ',
c     &              u_parallel, u_parallel_ratio, sign_u_parallel
c            if(id_num_FB_i.gt.0.and.id_num_FB_i.le.num_FB)then
c              print*,'+ F_fb_z_i ',F_fb_z_i
c              print*,' Z_Friction(id_num_FB_i) ', 
c     &                 Z_Friction(id_num_FB_i) 
c              print*,'nb_inFriction(id_num_FB_i) ', 
c     &                nb_inFriction(id_num_FB_i) 
c            endif
c            if(id_num_FB_j.gt.0.and.id_num_FB_j.le.num_FB)then
c              print*,'- F_fb_z_j ',-F_fb_z_j
c              print*,' Z_Friction(id_num_FB_j) ', 
c     &                 Z_Friction(id_num_FB_j) 
c              print*,'nb_inFriction(id_num_FB_j) ', 
c     &                nb_inFriction(id_num_FB_j) 
c            endif
c          endif
c        endif
c      endif

c         fxbp = 0.0    !For debugging
c         fzbp = 0.0    !For debugging
      end	