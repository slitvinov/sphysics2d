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


      subroutine monaghanBC(iip,jjp,drx,drz,dot,duxp,duzp,
     &                                    fxbp,fzbp)

      include 'common.2D'
      

      !print*,'monBC, iip,jjp ',iip,jjp

c     Boundary condition JWPCOE, 32
c
c     iip is a WaterParticle (WP)
c     jjp is a BoundaryParticle (BP)
c
c     drx = xp(WP) - xp(BedP)
c     drz = zp(WP) - zp(BedP)


c      Projection of WP position onto normal
       pdot3 = drx*xnb(jjp)+drz*znb(jjp)

c       if(pdot3.gt.0.0)then  !ensures that WP is on correct side of boundary
c
c         xpp = magnitude(vec(dr) x vec(n))   cross-product
c
c         xpp=(drx*ynb(jjp)-drz*xnb(jjp))


c        Projection of WP on plane of boundary
c
c        r(BP->WP) = [drx,drz,drz]
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
      
         !- Correction for particle exactly deltappt from another particle -
         g_dot_t = grx*xtb(jjp) + grz*ztb(jjp)
         xpp_orig = xpp
         if(.NOT.(abs(xpp).lt.deltappt))then  !i.e. >=deltappt
           if(g_dot_t.gt.0.0)then
             if(xpp.gt.0.0)then
               xpp = 0.999*xpp
             endif
           else   !if(g_dot_t.lt.0.0)then
             if(xpp.lt.0.0)then
               xpp = 0.999*xpp
             endif
           endif
         elseif(.NOT.(abs(xpp).gt.0.0))then  !i.e. ==0.0
           if(g_dot_t.gt.0.0)then
             xpp = xpp - 0.001*deltappt
           else   !if(g_dot_t.lt.0.0)then
             xpp = xpp + 0.001*deltappt
           endif
         endif
         !if(.NOT.(xpp.lt.deltappt))then
         !  xpp = xpp - 0.01*deltappt
         !end if  
         
         Prx = abs(xpp)/deltappt    
         
         facp = 0.0
         f_tb = 0.0
         f_fb = 0.0
         chi_X = 0.0
         Q_fn = 0.0
         i_CornerParticle = 0

         if(pdot3.gt.0.0)then ! particle inside
             q=pdot3 * one_over_2h           !normalized distance from wall
         else if(pdot3.lt.0.0)then !Particle ouside 
             q=0.2*abs(pdot3) * one_over_2h  !normalized distance from wall
         else !Particle on the boundary
             q=h*(1.0e-3) * one_over_2h      !correcting for zero distance from wall
         end if

c         if(q.lt.1.0.and.abs(xpp).lt.deltappt
c     &              .and.((iip.gt.nb.and.pdot3.gt.0.0).or.
c     &                    (iip.lt.nbp1)))then
         if(q.lt.1.0.and.((iip.gt.nb.and.pdot3.gt.0.0).or.
     &                    (iip.lt.nbp1)))then
     
             
c          u_normal < 0 if particles are approaching
c          u_normal > 0 if particles are moving apart
      
           !-- Corner Particle Identification --
           i_minus1 = iBP_Pointer_Info(iip,3)     !Absolute index of i-1 neighbour BP iip (Floating bodies only)
           i_plus1  = iBP_Pointer_Info(iip,4)     !Absolute index of i+1 neighbour BP iip (Floating bodies only)
           j_minus1 = iBP_Pointer_Info(jjp,3)     !Absolute index of i-1 neighbour BP jjp
           j_plus1  = iBP_Pointer_Info(jjp,4)     !Absolute index of i+1 neighbour BP jjp
           
           !- Calculation of normal velocity, Repulsive singularity function (Q_fn) -
           !- and tangential distance reduction function chi_X                      -
           if(xpp.lt.0.0.and.jjp.eq.j_minus1)then
             !Identify partner-CornerParticle index
             if((abs(xp(jjp-1)-xp(jjp)).lt.h_over_100).and.
     &          (abs(zp(jjp-1)-zp(jjp)).lt.h_over_100))then
               jjp_minus1 = jjp - 1
             elseif(jjp.gt.nbfm)then  !Floating Objects only
               id_num_FB_j = i_FB_Pointer_Info(jjp)
               if(id_num_FB_j.eq.1)then
                 i_FB_start = nbfm + 1
               else
                 i_FB_start = nb_FB(id_num_FB_j-1) + 1
               end if
               if((abs(xp(i_FB_start)-xp(jjp)).lt.h_over_100).and.
     &            (abs(zp(i_FB_start)-zp(jjp)).lt.h_over_100))then
                 jjp_minus1 = i_FB_start + 1
               else
                 jjp_minus1 = j_minus1
               endif
             else
               jjp_minus1 = j_minus1
             endif
             if(jjp_minus1.ne.j_minus1)then
               drp = sqrt(drx*drx + drz*drz)
               one_over_drp = 1.0/(drp + h_over_1000)
               q_r = drp*one_over_2h
               if(q_r.lt.1.0)then
                 Q_fn = (1.0-q_r)/sqrt(q_r)
               else
                 Q_fn = 0.0
               end if
               u_normal = (duxp*drx+duzp*drz)*one_over_drp
               chi_X = pdot3*one_over_drp
               vn1_dot_n2 =xnb(jjp)*xnb(jjp_minus1)
     &                   + znb(jjp)*znb(jjp_minus1)   !I need to know the indices of the other CP!!!
               vn1_dot_rij = (xnb(jjp)*drx + znb(jjp)*drz)*one_over_drp
               if(vn1_dot_rij.lt.vn1_dot_n2)then  !Activate only for WPs within arc of Corner Particles normals
                 Q_fn = 0.0
               end if
               xtb_r =  drz*one_over_drp
               ztb_r = -drx*one_over_drp
             else
               Q_fn = (1.0-q)/sqrt(q)
               u_normal = duxp*xnb(jjp) + duzp*znb(jjp)
               chi_X = (1.0 - Prx)
               q_r = q
               xtb_r = xtb(jjp)
               ztb_r = ztb(jjp)
             endif
             i_CornerParticle = 1
           else if(xpp.ge.0.0.and.jjp.eq.j_plus1)then
             !Identify partner-CornerParticle index
             if((abs(xp(jjp+1)-xp(jjp)).lt.h_over_100).and.
     &          (abs(zp(jjp+1)-zp(jjp)).lt.h_over_100))then
               jjp_plus1 = jjp + 1
             elseif(jjp.gt.nbfm)then  !Floating Objects only
               id_num_FB_j = i_FB_Pointer_Info(jjp)
               if(id_num_FB_j.eq.1)then
                 i_FB_start = nbfm + 1
               else
                 i_FB_start = nb_FB(id_num_FB_j-1) + 1
               end if
               if((abs(xp(i_FB_start)-xp(jjp)).lt.h_over_100).and.
     &            (abs(zp(i_FB_start)-zp(jjp)).lt.h_over_100))then
                 jjp_plus1 = i_FB_start + 1
               else
                 jjp_plus1 = j_plus1
               endif
             else
               jjp_plus1 = j_plus1
             endif
             if(jjp_plus1.ne.j_plus1)then
               drp = sqrt(drx*drx + drz*drz)
               one_over_drp = 1.0/(drp + h_over_1000)
               q_r = drp*one_over_2h
               if(q_r.lt.1.0)then
                 Q_fn = (1.0-q_r)/sqrt(q_r)
               else
                 Q_fn = 0.0
               end if
               u_normal = (duxp*drx+duzp*drz)*one_over_drp
               chi_X = pdot3*one_over_drp
               
               vn1_dot_n2 = xnb(jjp)*xnb(jjp_plus1)
     &                    + znb(jjp)*znb(jjp_plus1)   !I need to know the indices of the other CP!!!
               vn1_dot_rij = (xnb(jjp)*drx + znb(jjp)*drz)*one_over_drp
               if(vn1_dot_rij.lt.vn1_dot_n2)then  !Activate only for WPs within arc of Corner Particles normals
                 Q_fn = 0.0
               end if
               xtb_r =  drz*one_over_drp
               ztb_r = -drx*one_over_drp
             else
               Q_fn = (1.0-q)/sqrt(q)
               u_normal = duxp*xnb(jjp) + duzp*znb(jjp)
               chi_X = (1.0 - Prx)
               q_r = q
               xtb_r = xtb(jjp)
               ztb_r = ztb(jjp)
             endif
             i_CornerParticle = 1
           else if(abs(Prx).lt.1.0)then  !Conventional Boundary Particle
             Q_fn = (1.0-q)/sqrt(q)
             u_normal = duxp*xnb(jjp) + duzp*znb(jjp)
             chi_X = (1.0 - Prx)
             q_r = q
             xtb_r = xtb(jjp)
             ztb_r = ztb(jjp)
           end if
             
             
           if(chi_X.gt.0.0.and.Q_fn.gt.0.0)then
             !c0 = coeff.Vmax (Vmax = sqrt(g.h_SWL); coef = beta_coef.10, beta_coef2 = beta_coef**2
             !beta_coef2 = Factor that B is changed by 
             if(iip.gt.nb)then  !Normal FP-BP Interaction
               cs_i = cs0  !cs(iip)
               if(u_normal.gt.0.0)then
                 App=(0.01*cs_i*cs_i*one_over_beta_coef2
     &                 -cs_i*u_normal*one_over_beta_coef)*one_over_h
               else
                 App=0.01*(cs_i*cs_i*one_over_beta_coef2)
     &                   *one_over_h    
               endif
             else !BP-BP interactions
               cs_i = cs0
               if(u_normal.gt.0.0)then    !Particles receding:  u_normal > 0
                 !cs2_diff = 0.5*abs(cs_i*cs_i)   !
                 cs2_diff = 1.0*abs(cs_i*cs_i)   !
                 App=(0.01*cs2_diff*one_over_beta_coef2
     &               -cs_i*u_normal*one_over_beta_coef)*one_over_h
                 q2 = q_r*q_r
                 f_of_q = (-2.0*(q_r*q2) + 3.0*(q2))
                 repulsion_factor = 10.0
     &                          + (10.-10.)*(1.0-f_of_q)/sqrt(f_of_q)
                 App=(0.01*cs2_diff*one_over_beta_coef2
     &                 +repulsion_factor*cs_i*u_normal
     &                 *one_over_beta_coef)*one_over_h
               else                       !Particles approaching:  u_normal < 0
                 cs2_diff = abs(cs_i*cs_i)   !Set this for elastic bouncing
                 !cs2_diff = 0.0   !Set this for inelastic bouncing
                 !repulsion_factor = 10.0
                 q2 = q_r*q_r
                 f_of_q = (-2.0*(q_r*q2) + 3.0*(q2))
                 repulsion_factor = 10.0
     &                          + (1000.-10.)*(1.0-f_of_q)/sqrt(f_of_q)
                 App=(0.01*cs2_diff*one_over_beta_coef2
     &                 -repulsion_factor*cs_i*u_normal
     &                 *one_over_beta_coef)*one_over_h
               endif
               App = App*0.5
               !if(q.lt.0.2.and.(jjp.ne.i_minus1.and.jjp.ne.i_plus1))then
                 q = q_r*q_r
               !end if
               if(iip.eq.i_plus1.or.iip.eq.i_minus1)then
                 App = App*0.5  !Halves contribution for corner particles placed ontop of each other
               endif
             endif

             !if(iip.lt.nbp1)epsilon_Z = 1.0  !Has no effect on the block's Vertical position
       
c          --- Depth Correction ---
             eps_min = 0.02
             epsilon_Z = abs(1.0 - (zp(jjp)/h_SWL))
             dyn_head =(up(jjp)*up(jjp) + wp(jjp)*wp(jjp))
     &                  /(2.0*abs(grz)*h_SWL)
c     &                  /(2.0*9.81*h_SWL)   !1-D
           
             if(epsilon_Z.gt.1.0-eps_min)then !Close to the bottom (No correction )
                epsilon_Z = 1.0
             elseif(zp(iip).gt.h_SWL)then !over still water level          
                epsilon_Z = eps_min + dyn_head
             else                 !under still water level & far from bottom
                epsilon_Z = eps_min + epsilon_Z + dyn_head
             endif

c          --- Velocity Correction ---
             if(u_normal.gt.0.0)then
                epsilon_dyn = 0.0
             else 
                epsilon_dyn = min(1.0,40.*abs(u_normal)/cs0 )    !factor 20 needed 
                                                       ! to keep scheme stable                                                                                  ! with Monaghan BC wavemaker
             endif

             !if(epsilon_dyn.gt.1.0)then
             !   epsilon_dyn = 1.0
             !endif
             
             if(iip.lt.nbp1)then
                epsilon_Z = 1.0
             endif

             Appc = App*(epsilon_Z + epsilon_dyn) !Force correction
             
c             facp=0.5*(1.+cos(pi*abs(xpp)/deltappt))*
c     1            (Appc*0.5*(1.0-q)/sqrt(q))
c     2            *vnorm_mass

             !Prx = abs(xpp)/deltappt
             !chi_X = (1.0 - Prx)
             !Q_fn = (1.0-q)/sqrt(q)
             fmax = Appc*Q_fn*vnorm_mass
             facp = chi_X*fmax
             
             fxbp = xnb(jjp)*facp
             fzbp = znb(jjp)*facp   

             if(iip.gt.nb.and.jjp.lt.nbp1)then  !Floating Body Version
c              !-- Angular momentum correction --
               g_dot_t = grx*xtb_r + grz*ztb_r  !Be careful of this:  g_dot_t = grx*xtb(jjp) + grz*ztb(jjp)
               g_dot_t_normalised = abs(g_dot_t)/abs(grz)
               if(xpp.gt.0.0)then
                 f_tb =  two_alpha_h*g_dot_t_normalised*fmax*
     &                   (0.25*deltappt*one_over_r_ij_dot_n)
               else
                 f_tb = -two_alpha_h*g_dot_t_normalised*fmax*
     &                   (0.25*deltappt*one_over_r_ij_dot_n)
               end if
               
               !- Correction for WP near bottom corner or surface
               h_SWL_minus_zp = h_SWL - zp(iip)
               !rNum_h = 3.0*h
               !one_over_rNum_h = 1.0/rNum_h
               f_tb_orig = f_tb
               if(h_SWL_minus_zp.lt.0.0)then  !WP is above h_SWL
                 f_tb = 0.0  !*f_tb
               elseif(h_SWL_minus_zp.lt.rNum_h)then
                 epsilon_Z_AMC = h_SWL_minus_zp * one_over_rNum_h
                 if(epsilon_Z_AMC.lt.0.90)then
                   epsilon_Z_AMC = (epsilon_Z_AMC + 0.10)
                 endif
                 f_tb = f_tb * epsilon_Z_AMC
               !elseif(zp(iip).lt.rNum_h)then  !i.e. less than 2h above the bottom
               !  epsilon_Z_AMC = zp(iip) * one_over_rNum_h
               !  f_tb = f_tb * epsilon_Z_AMC
               endif

               !fxbp = fxbp + xtb(jjp)*f_tb
               !fzbp = fzbp + ztb(jjp)*f_tb
               fxbp = fxbp + xtb_r*f_tb
               fzbp = fzbp + ztb_r*f_tb
             endif
             !abs_t_dot_t = abs(xtb(iip)*xtb(jjp) + ztb(iip)*ztb(jjp))
             abs_t_dot_t = abs(xtb(iip)*xtb_r + ztb(iip)*ztb_r)
             if(iip.lt.nbp1)then
               !fmax_t = App_t*Q_fn*vnorm_mass
               !facp_t = chi_X*fmax_t
               !-- Friction force between rigid objects --
               id_num_FB_i = i_FB_Pointer_Info(iip)
               pmu = friction_coeff(id_num_FB_i)
               if(abs_t_dot_t.lt.0.6)then  !Activates forces only between sliding faces
                 pmu = 0.0
                 fxbp = 0.0
                 fzbp = 0.0

                 facp_F = 0.0
               else
                 facp_F = fmax  !facp  !_t
               endif
               if(jjp.gt.nbfm)pmu = 0.0  !No friction between moving blocks   !<--------------------
               !u_parallel = duxp*xtb(jjp) + duzp*ztb(jjp)
               u_parallel = duxp*xtb_r + duzp*ztb_r
               !f_tb = - sign(1.0,u_parallel)*pmu*facp
               if(u_parallel.gt.0.0)then
                 sign_u_parallel = 1.0
               elseif(u_parallel.lt.0.0)then
                 sign_u_parallel = -1.0
               else
                 sign_u_parallel = 0.0
               endif
               !f_fb = - sign_u_parallel*pmu*facp
               
c               f_fb = - sign_u_parallel*pmu*facp_F
c               F_Dynamic = facp_F     
c               f_fb_min = chi_X*Appc*vnorm_mass  !(i.e. same as F_Dynamic but without singularity Q_fn)               
c               f_fb = - sign_u_parallel*pmu*min(f_fb_min, F_Dynamic)   

               if(num_FB.lt.2)then
                 !F_Static = abs(X_nonFriction(1))
                 F_Static = (X_nonFriction(1)*xtb_r
     &                     + Z_nonFriction(1)*ztb_r)
                 F_Static_abs = abs(F_Static)
               else
                 F_Static = (X_nonFriction(1)*xtb_r
     &                     + Z_nonFriction(1)*ztb_r)
                 F_Static_abs = abs(F_Static)
               endif
               if(F_Static.gt.0.0)then
                 sign_F_Static = 1.0
               elseif(F_Static.lt.0.0)then
                 sign_F_Static = -1.0
               else
                 sign_F_Static = 0.0
               endif
               if(sign_u_parallel.eq.sign_F_Static)then
                 !F_Static_abs = F_Static_abs * 1.0005   !An extra 2.5% for cases where the block inches forward
               endif
               
c               fxbp = fxbp + f_fb*xtb_r   !+ f_fb*xtb(jjp)
c               fzbp = fzbp + f_fb*ztb_r   !+ f_fb*ztb(jjp)
               
               if(iip.gt.nbfm)then
                 id_num_FB_i = i_FB_Pointer_Info(iip)
                 F_Dynamic = pmu*facp_F*pm(iip)
                 u_parallel_ratio = abs(u_parallel)
     &                             *one_over_u_parallel_max
                 if(u_parallel_ratio.lt.1.0)then
                   !- Static Friction Force -
                   f_temp = - sign_u_parallel*chi_X
     &             *min(F_Dynamic,F_Static_abs)  !ensures max friction force is <= F_Dynamics
c     &             *min(F_Dynamic*u_parallel_ratio,F_Static_abs)
                   F_fb_x =   f_temp*xtb_r
                   F_fb_z =   f_temp*ztb_r
                 elseif(u_parallel_ratio.lt.2.0)then
                   !- Static-to-Dynamic Friction Force F_SD -
                   F_SD = F_Static_abs + (F_Dynamic - F_Static_abs)
     &                           *(u_parallel_ratio - 1.0)  !linear interpolation
                   f_temp = - sign_u_parallel*chi_X
     &             *min(F_Dynamic,F_SD)  !ensures max friction force is <= F_Dynamics
c     &             *min(F_Dynamic*u_parallel_ratio,F_Static_abs)
                   F_fb_x =   f_temp*xtb_r
                   F_fb_z =   f_temp*ztb_r
                 else
                   !- Dynamic Friction Force -
                   f_temp = - sign_u_parallel*chi_X*F_Dynamic
                   F_fb_x =   f_temp*xtb_r
                   F_fb_z =   f_temp*ztb_r
                 endif
                 X_Friction(id_num_FB_i) = X_Friction(id_num_FB_i)
     &                                        + F_fb_x
                 Z_Friction(id_num_FB_i) = Z_Friction(id_num_FB_i)
     &                                        + F_fb_z
                 fxbp = fxbp + F_fb_x/pm(iip)  !per unit mass force
                 fzbp = fzbp + F_fb_z/pm(iip)  !per unit mass force
                 if(chi_X.gt.0.5.and.abs_t_dot_t.gt.0.6)then
                   nb_inFriction(id_num_FB_i) = 
     &             nb_inFriction(id_num_FB_i) + 1
                 endif
               endif
               !- This section needs to be updated for multiple floating bodies sliding past each other -
               if(jjp.gt.nbfm.and.iip.lt.nbp1)then
                 id_num_FB_j = i_FB_Pointer_Info(jjp)
                 F_Dynamic = pmu*facp_F*pm(iip)
                 F_fb_x = - sign_u_parallel*chi_X*xtb_r
     &                     *min(F_Dynamic,F_Static_abs)
                 F_fb_z = - sign_u_parallel*chi_X*ztb_r
     &                     *min(F_Dynamic,F_Static_abs)
                 X_Friction(id_num_FB_j) = X_Friction(id_num_FB_j)
     &                                        - F_fb_x      ! - f_fb*xtb_r*pm(jjp)  !CHECK MINUS SIGN HERE!!!!!
                 bigWdot(id_num_FB_j) = bigWdot(id_num_FB_j)
     &                                        - F_fb_z      !- f_fb*ztb_r*pm(jjp)
                 if(chi_X.gt.0.5.and.abs_t_dot_t.gt.0.6)then
                   nb_inFriction(id_num_FB_j) = 
     &             nb_inFriction(id_num_FB_j) + 1
c                   print*,'i,j,xp,zp, nb_inFriction#',
c     &             iip,jjp,xp(iip),zp(iip),xp(jjp),zp(jjp),
c     &             nb_inFriction(1)
                 endif
               endif
                          
             end if    !End of:  if(1.eq.2.and.iip.gt.nb)then
           else  
             fxbp = 0.0
             fzbp = 0.0
           endif   !End of:  if(chi_X.gt.0.0.and.Q_fn.gt.0.0)then
           
         else
           fxbp = 0.0
           fzbp = 0.0
         end if
      end     
