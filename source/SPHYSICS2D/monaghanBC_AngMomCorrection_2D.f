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

         chi_X = 0.0
         Q_fn = 0.0
         
         if(xpp.gt.0.0)then
           deltappt = deltaptb(jjp,2)
         else
           deltappt = deltaptb(jjp,1)
         end if
         
         !- Correction for particle exactly deltappt from another particle -
         g_dot_t = grx*xtb(jjp) + grz*ztb(jjp)
         if(.NOT.(abs(xpp).lt.deltappt))then  !i.e. >=deltappt
           if(g_dot_t.gt.0.0)then
             if(xpp.gt.0.0)then
               xpp = 0.99*xpp
             endif
           else   !if(g_dot_t.lt.0.0)then
             if(xpp.lt.0.0)then
               xpp = 0.99*xpp
             endif
           endif
         elseif(.NOT.(abs(xpp).gt.0.0))then  !i.e. ==0.0
           if(g_dot_t.gt.0.0)then
             xpp = xpp - 0.01*deltappt
           else   !if(g_dot_t.lt.0.0)then
             xpp = xpp + 0.01*deltappt
           endif
         endif
         !if(.NOT.(xpp.lt.deltappt))then
         !  xpp = xpp - 0.01*deltappt
         !end if      

         if(pdot3.gt.0.0)then ! particle inside
             q=pdot3 * one_over_2h           !normalized distance from wall
         else if(pdot3.lt.0.0)then !Particle ouside 
             q=0.2*abs(pdot3) * one_over_2h  !normalized distance from wall
         else !Particle on the boundary
             q=h*(1.0e-3) * one_over_2h      !correcting for zero distance from wall
         end if

         if(q.lt.1.0.and.pdot3.gt.0.0.and.
     &                abs(xpp).lt.deltappt)then
     
             u_normal = duxp*xnb(jjp) + duzp*znb(jjp)
c            u_normal < 0 if particles are approaching
c            u_normal > 0 if particles are moving apart
      
             !c0 = coeff.Vmax (Vmax = sqrt(g.h_SWL); coef = beta_coef.10, beta_coef2 = beta_coef**2
             !beta_coef2 = Factor that B is changed by 
             if(u_normal.gt.0.0)then
              App=(0.01*cs(iip)*cs(iip)*one_over_beta_coef2
     &                -cs(iip)*u_normal*one_over_beta_coef)*one_over_h
             else
              App=0.01*(cs(iip)*cs(iip)*one_over_beta_coef2)*one_over_h    
             endif

       
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
                epsilon_dyn = 40.*abs(u_normal)/cs0     !factor 20 needed 
                                                       ! to keep scheme stable                                                                                  ! with Monaghan BC wavemaker
             endif

             if(epsilon_dyn.gt.1.0)then
                epsilon_dyn = 1.0
             endif

             Appc = App*(epsilon_Z + epsilon_dyn) !Force correction
             
c             facp=0.5*(1.+cos(pi*abs(xpp)/deltappt))*
c     1            (Appc*0.5*(1.0-q)/sqrt(q))
c     2            *vnorm_mass

             Prx = abs(xpp)/deltappt
             chi_X = (1.0 - Prx)
             Q_fn = (1.0-q)/sqrt(q)
             fmax = Appc*Q_fn  !*vnorm_mass
             facp = chi_X*fmax
             
             fxbp = xnb(jjp)*facp
             fzbp = znb(jjp)*facp   
         
             if(1.eq.1)then
             !if(1.eq.2.and.iip.gt.nb.and.jjp.lt.nbp1)then  Floating Body Version
               !-- Angular momentum correction --

c               if(abs(pdot3).lt.eta)then   !Singularity catching
c                 r_ij_dot_n = eta
c               else if(pdot3.lt.0.0)then   !Particle wrong side of boundary
c                 r_ij_dot_n = abs(pdot3)
c               else                        !Normal behaviour
c                 r_ij_dot_n = pdot3
c               end if
c               r_ij_dot_n = h/1.3   !dpx
c
c               two_alpha_h = 2.0*1.3
c               g_dot_t = abs(grx*xtb(jjp) + grz*ztb(jjp))/abs(grz)
c               if(xpp.gt.0.0)then
c                 f_tb =  two_alpha_h*g_dot_t*fmax*
c     &                   (0.25*deltappt/r_ij_dot_n)
c               else
c                 f_tb = -two_alpha_h*g_dot_t*fmax*
c     &                   (0.25*deltappt/r_ij_dot_n)
c               end if
c               
c               !- Correction for WP near bottom corner or surface
c               h_SWL_minus_zp = h_SWL - zp(iip)
c               rNum_h = 3.0*h
c               one_over_rNum_h = 1.0/rNum_h
c               !f_tb_orig = f_tb
c               if(h_SWL_minus_zp.lt.0.0)then  !WP is above h_SWL
c                 f_tb = 0.10*f_tb
c               elseif(h_SWL_minus_zp.lt.rNum_h)then
c                 epsilon_Z_AMC = h_SWL_minus_zp * one_over_rNum_h
c                 if(epsilon_Z_AMC.lt.0.90)then
c                   epsilon_Z_AMC = (epsilon_Z_AMC + 0.10)
c                 endif
c                 f_tb = f_tb * epsilon_Z_AMC
c               !elseif(zp(iip).lt.rNum_h)then  !i.e. less than 2h above the bottom
c               !  epsilon_Z_AMC = zp(iip) * one_over_rNum_h
c               !  f_tb = f_tb * epsilon_Z_AMC
c               endif
               
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
               
               fxbp = fxbp + xtb(jjp)*f_tb
               fzbp = fzbp + ztb(jjp)*f_tb
           
             endif
         
         else
             fxbp = 0.0
             fzbp = 0.0
         end if

  
      end     
