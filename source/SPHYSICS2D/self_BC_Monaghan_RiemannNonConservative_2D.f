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

      subroutine self(j1,kind_p1,ini_kind_p2)
c
      include 'common.2D'
      

      do kind_p2=ini_kind_p2,2
    
      jj_start = 1

      if(kind_p1.gt.kind_p2)then
        jj_start = nplink_max + 1      
      endif
      
      do ii=1,nc(j1,kind_p1)
         i = ibox(j1,kind_p1,ii)
        
         if(kind_p1.eq.kind_p2)then
           jj_start = jj_start + 1
         end if
         
         do jj=jj_start,nc(j1,kind_p2)
            j = ibox(j1,kind_p2,jj)

           drx = xp(i) - xp(j)
           drz = zp(i) - zp(j)

           rr2 = drx*drx + drz*drz

           if(rr2.lt.fourh2.and.rr2.gt.1.e-18) then

             dux = up(i) - up(j)
             duz = wp(i) - wp(j)


c            Calculating kernel & Normalized Kernel Gradient
             call kernel(drx,drz,i,j,j1,j1,rr2) 
             call kernel_correction(i,j)

c
c            ...  average density
c
             robar  = 0.5*( rhop(i) + rhop(j) )
             one_over_rhobar = 2.0/(rhop(i) + rhop(j))
c
             cbar = 0.5*( cs(i)   + cs(j)   ) 
c
c            ...  inner product rv
c
             dot = drx*dux + drz*duz

c            Used to calculate the time step due to viscosity

             visc_dt=max(dot/(rr2 + eta2),visc_dt) 
c
c            ...  pressure and viscous force (1992; 3.3)
c
c            pm(j) is mass of particle j
c
c             p_v = pr(i) + pr(j)


             if(index_tensile.eq.1) then   

c              _____ Tensile correction

              fab=Wab*od_Wdeltap                 !NOTE: We'll use a non-normalized
              fab=fab*fab                        !kernel to calculate tensile correction
              fab=fab*fab                        !  It's the (Wab/Wdeltap)**4  of Monaghan's paper

              if (p(i).gt.0) then
                  Ra= 0.01 *pr(i)
              else
                  Ra= 0.2 *abs(pr(i))
              endif

              if (p(j).gt.0) then
                  Rb= 0.01 *pr(j)
              else
                  Rb= 0.2 *abs(pr(j))
              endif

              !p_v = p_v+ (Ra+Rb)*fab

            endif

		if(itime.ge.29.and.i.eq.8496.and.j.eq.8504)then
	write(146,*) 'antes',itime,rhop(i), rhop(j)
	endif
            if(i.gt.nb.and.j.gt.nb) then   !both fluid particles

              call approx_RiemannSolver(i,j,drx,drz,rr2,cbar,robar,
     &                              u_a,w_a,rho_a,p_a,TE_a)
              frx_tmpj = pVol(j)*2.0*frxi
              frz_tmpj = pVol(j)*2.0*frzi

              frx_tmpi = pVol(i)*2.0*frxj
              frz_tmpi = pVol(i)*2.0*frzj

              p_v_x1 = p_a
              p_v_z1 = p_a

              ax(i) = ax(i) - p_v_x1 * frx_tmpj * one_over_rhop(i)
              az(i) = az(i) - p_v_z1 * frz_tmpj * one_over_rhop(i)

              ax(j) = ax(j) + p_v_x1 * frx_tmpi * one_over_rhop(j)
              az(j) = az(j) + p_v_z1 * frz_tmpi * one_over_rhop(j)
              

              call gradients_calc(i,j,dux,duz)               

              call viscosity(dot,drx,drz,dux,duz,rr2,
     +              cbar,robar,one_over_rhobar,i,j,j1,j1,term2i,term2j)

c             ...  density acceleration
c             using the derivative of the kernel, not the kernel itself
              
              dot2ai = (u_a - up(j))*frx_tmpi
     &               + (w_a - wp(j))*frz_tmpi
              dot2aj = (u_a - up(i))*frx_tmpj
     &               + (w_a - wp(i))*frz_tmpj
              ar(i) = ar(i) - rhop(i)*dot2aj
              ar(j) = ar(j) + rhop(j)*dot2ai

c             ...  Thermal Energy

              aTE(i) = aTE(i) - p_a*dot2aj * one_over_rhop(i)
              aTE(j) = aTE(j) + p_a*dot2ai * one_over_rhop(j)

c
c             ...  XSPH correction (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.6)
c
              pmj_Wab_over_rhobar = pm(j)*Wab*one_over_rhobar
              ux(i) = ux(i) - dux*pmj_Wab_over_rhobar  !pm(j) * dux * Wab / robar
              wx(i) = wx(i) - duz*pmj_Wab_over_rhobar  !pm(j) * duz * Wab / robar

        
              pmi_Wab_over_rhobar = pm(i)*Wab*one_over_rhobar
              ux(j) = ux(j) + dux*pmi_Wab_over_rhobar   !pm(i) * dux * Wab / robar
              wx(j) = wx(j) + duz*pmi_Wab_over_rhobar   !pm(i) * duz * Wab / robar


c             ...  Vorticity calculation

            if(ipoute.eq.1.and.i_vort.eq.1.and.i.gt.nb.and.j.gt.nb)then
                  call vorticity_calc(i,j,dux,duz)
            endif 

           else if(i.gt.nb.and.j.le.nb)then              !j is boundary particle
             call MonaghanBC(i,j,drx,drz,dot,dux,duz,
     &                       fxbMON,fzbMON)    
             ax(i) = ax(i) + fxbMON
             az(i) = az(i) + fzbMON

              !- Viscous wall force -
             temp = pVol(i)*visc_wall*
     &         ((drx*frxi+drz*frzi)/(rr2 + eta2))
             ax(i) = ax(i) + temp*dux
             az(i) = az(i) + temp*duz


           else if(j.gt.nb.and.i.le.nb)then                  !i is boundary particle
             call MonaghanBC(j,i,-drx,-drz,dot,-dux,-duz,
     &                       fxbMON,fzbMON)   
             ax(j) = ax(j) + fxbMON
             az(j) = az(j) + fzbMON

             !- Viscous wall force -
             temp = pVol(j)*visc_wall*
     &         ((drx*frxj+drz*frzj)/(rr2 + eta2))
             ax(j) = ax(j) - temp*dux
             az(j) = az(j) - temp*duz


           endif ! Interaction fluid-fluid or fluid-boundary
         endif ! if q<2
        enddo ! Box jj
       enddo   ! Box ii
      enddo   ! Kind of particle

      end
c
