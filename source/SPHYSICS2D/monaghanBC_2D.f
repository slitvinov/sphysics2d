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
       pdot3 = drx*xnb(jjp)+drz*znb(jjp)      !SPHysics Guide, psi in Eq. 2.13

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

c        Projection of plane-projected vector r(BP->WP projection) onto tangent     !SPHysics Guide, ksi in Eq. 2.12 & 2.15

         xpp = drx_proj*xtb(jjp)+drz_proj*ztb(jjp)

         if(xpp.gt.0.0)then
           deltappt = deltaptb(jjp,2)     !SPHysics Guide, delta_b in Eq. 2.15
         else
           deltappt = deltaptb(jjp,1)     !SPHysics Guide, delta_b in Eq. 2.15
         end if
      

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
              App=0.01*(cs(iip)*cs(iip)*one_over_beta_coef2)*one_over_h     !SPHysics Guide, Eq. 2.14  
             endif

	
c          --- Depth Correction ---      !SPHysics Guide, Eq. 2.17
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

c          --- Velocity Correction ---     !SPHysics Guide, Eq. 2.18
             if(u_normal.gt.0.0)then
                epsilon_dyn = 0.0
             else 
                epsilon_dyn = 40.*abs(u_normal)/cs0     !factor 20 needed  
			                                  ! to keep scheme stable 											      ! with Monaghan BC wavemaker
             endif

             if(epsilon_dyn.gt.1.0)then
                epsilon_dyn = 1.0
             endif

             Appc = App*(epsilon_Z + epsilon_dyn) !Force correction     !SPHysics Guide, Eq. 2.16
             
             facp=0.5*(1.+cos(pi*abs(xpp)/deltappt))*
     1            (Appc*(1.0-q)/sqrt(q))       !SPHysics Guide, Eq. 2.12, 2.13, 2.15
c     2            *vnorm_mass

             fxbp = xnb(jjp)*facp   !SPHysics Guide, Eq. 2.12
             fzbp = znb(jjp)*facp   !SPHysics Guide, Eq. 2.12

         else
             fxbp = 0.0
             fzbp = 0.0
         end if

        

      end	
