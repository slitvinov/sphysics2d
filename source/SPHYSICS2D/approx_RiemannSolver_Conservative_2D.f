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

      subroutine approx_RiemannSolver(i,j,drx,drz,rr2,
     &                         drho,drhou,drhow,dTotalE,
     &                         u_a,w_a,rho_a,p_a,E_a,u_0,w_0)
c
      include 'common.2D'
      
      real(kind=kind_XZ ) x_ij_av, z_ij_av
      real(kind=kind_XZ ) dx_i, dz_i, dx_j, dz_j

c     --- HLLC Riemann Solver - see Toro (1999) ---

      one_over_rad = 1.0/sqrt(rr2)

      sin_theta = -drz*one_over_rad
      cos_theta = -drx*one_over_rad
           
      if(iTVD.ne.1)then  !First-order scheme
        u_ij = up(i)
        w_ij = wp(i)

        u_ji = up(j)
        w_ji = wp(j)

        rho_ij = rhop(i)
        rho_ji = rhop(j)
        E_ij = TotalE(i)
        E_ji = TotalE(j)
      else               !Higher-Order using MUSCL Reconstruction
        x_ij_av = 0.5*(xp(i) + xp(j))
c       y_ij_av = 0.5*(yp(i) + yp(j))
        z_ij_av = 0.5*(zp(i) + zp(j))
        dx_i = x_ij_av - xp(i)
c       dy_i = y_ij_av - yp(i)
        dz_i = z_ij_av - zp(i)
        dx_j = -(x_ij_av - xp(j))
c       dy_j = -(y_ij_av - yp(j))
        dz_j = -(z_ij_av - zp(j))
        delrhou_i = drhoudx_CSPHo(i)*dx_i
     &            + drhoudz_CSPHo(i)*dz_i
c       delrhov_i = drhovdx_CSPHo(i)*dx_i
c    &            + drhovdz_CSPHo(i)*dz_i
        delrhow_i = drhowdx_CSPHo(i)*dx_i
     &            + drhowdz_CSPHo(i)*dz_i
        delrho_i  = drhodx_CSPHo(i)*dx_i
     &            + drhodz_CSPHo(i)*dz_i
        delrhou_j = drhoudx_CSPHo(j)*dx_j
     &            + drhoudz_CSPHo(j)*dz_j
c       delrhov_j = drhovdx_CSPHo(j)*dx_j
c    &            + drhovdz_CSPHo(j)*dz_j
        delrhow_j = drhowdx_CSPHo(j)*dx_j
     &            + drhowdz_CSPHo(j)*dz_j
        delrho_j  = drhodx_CSPHo(j)*dx_j
     &            + drhodz_CSPHo(j)*dz_j
        delE_i = dTotalEdx_CSPHo(i)*dx_i
     &         + dTotalEdz_CSPHo(i)*dz_i
        delE_j = dTotalEdx_CSPHo(j)*dx_j
     &         + dTotalEdz_CSPHo(j)*dz_j
        d_rho_ij     = -drho      !rhop(j) - rhop(i)
        d_rhou_ij    = -drhou     !Volrhou(j)/pVol(j) - Volrhou(i)/pVol(i)
        d_rhow_ij    = -drhow     !Volrhow(j)/pVol(j) - Volrhow(i)/pVol(i)

        d_E_ij  = -dTotalE   !TotalE(j) - TotalE(i)
        !- MUSCL Limiters - Toro (2001) p. 208
        call limiter(i,j,delrhou_i,delrhou_j,d_rhou_ij,
     &                   del_rhou_lim_i,del_rhou_lim_j)
        call limiter(i,j,delrhow_i,delrhow_j,d_rhow_ij,
     &                   del_rhow_lim_i,del_rhow_lim_j)
        call limiter(i,j,delE_i,delE_j,d_E_ij,
     &                   del_E_lim_i,del_E_lim_j)
        call limiter(i,j,delrho_i,delrho_j,d_rho_ij,
     &                   del_rho_lim_i,del_rho_lim_j)

        rho_ij = rhop(i) + del_rho_lim_i
        rho_ji = rhop(j) - del_rho_lim_j
        u_ij = up(i) + del_rhou_lim_i/rho_ij
        w_ij = wp(i) + del_rhow_lim_i/rho_ij

        u_ji = up(j) - del_rhou_lim_j/rho_ji
        w_ji = wp(j) - del_rhow_lim_j/rho_ji

        E_ij = TotalE(i) + del_E_lim_i
        E_ji = TotalE(j) - del_E_lim_j
      end if

      u_L =  u_ij*cos_theta + w_ij*sin_theta
      w_L = -u_ij*sin_theta + w_ij*cos_theta

      u_R =  u_ji*cos_theta + w_ji*sin_theta
      w_R = -u_ji*sin_theta + w_ji*cos_theta

      rho_L = rho_ij
      rho_R = rho_ji
      E_L = E_ij
      E_R = E_ji
      u_0 = 0.5*(u_ij + u_ji)
      w_0 = 0.5*(w_ij + w_ji)

      u_0_R = u_0*cos_theta + w_0*sin_theta

c     --- Initial estimates of pressure, velocity and density in star region ---      
      
      if(i_EoS.ne.2)then  !Incompressible water
        call equation_of_state(rho_L,E_L,p_L,cs_L)
        call equation_of_state(rho_R,E_R,p_R,cs_R)
        if(i_EoS.eq.1)then
          !- Tait EoS -
          u_star  = 0.5*(u_L  + u_R ) + (cs_L - cs_R)/(gamma - 1.)
          cs_star = 0.5*(cs_L + cs_R) + 0.25*(gamma - 1.)*(u_L - u_R)
          rho_star_over_rho0 = (cs_star/cs0)**(1.0/3.0)
          rho_star = rho_star_over_rho0 * rho0
          p_star = B * ((rho_star_over_rho0)**i_gamma - 1.)
        elseif(i_EoS.eq.3)then
          !- Morris EoS -
          u_star = (u_R*rho_R*cs_R + u_L*rho_L*cs_L - p_R + p_L)
     &                  /(rho_L*cs_L + rho_R*cs_R)
          p_star = (p_R*rho_L*cs_L + p_L*rho_R*cs_R - 
     &                   rho_L*cs_L*rho_R*cs_R*(u_R - u_L))
     &                  /(rho_L*cs_L + rho_R*cs_R)
          rho_star = ((p_star - 0.0)/(cs0**2)) + rho0
          p_PVRS = p_star
        endif
      else  !Ideal Gas equation of state
        TE_L = E_L/rho_L - 0.5*(u_L*u_L + w_L*w_L)  !Left  Thermal Energy
        TE_R = E_R/rho_R - 0.5*(u_R*u_R + w_R*w_R)  !Right Thermal Energy
        call equation_of_state(rho_L,TE_L,p_L,cs_L)
        call equation_of_state(rho_R,TE_R,p_R,cs_R)
        c_hat = 0.5*(cs_L + cs_R)
        rho_hat = 0.5*(rho_L + rho_R)
        p_star = 0.5*(p_L + p_R) + 0.5*(u_L - u_R)*c_hat*rho_hat
c        u_star = 0.5*(u_L + u_R) + 0.5*(p_L - p_R)/(c_hat*rho_hat)
        !u_star  = 0.5*(u_L  + u_R ) + (cs_L - cs_R)/(gamma - 1.)
         p_star = (p_R*rho_L*cs_L + p_L*rho_R*cs_R - 
     &                 rho_L*cs_L*rho_R*cs_R*(u_R - u_L))
     &                /(rho_L*cs_L + rho_R*cs_R)
         u_star = (u_R*rho_R*cs_R + u_L*rho_L*cs_L - p_R + p_L)
     &                /(rho_L*cs_L + rho_R*cs_R)
c        if(1.eq.2)then
c         !Toro (1999) Section 10.5.2 & 9.5.2 for determining p* and u*
c         Q_user = 2.0
c         p_max = max(p_L,p_R)
c         p_min = min(p_L,p_R)
c         p_PVRS = 0.5*(p_L + p_R) + 0.5*(u_L - u_R)*c_hat*rho_hat
c         Q_pressureRatio = p_max/p_min
c         if(Q_pressureRatio.lt.Q_user.and.
c     &     (p_min.lt.p_PVRS.and.p_PVRS.lt.p_max))then
c           p_star = p_PVRS
c           u_star = 0.5*(u_L + u_R) + 0.5*(p_L - p_R)/(c_hat*rho_hat)         
c         else
c           if(p_PVRS.le.p_min)then  !Use Two-Rarefaction Riemann Solver (TRRS)
c             gamma_1 = (gamma - 1.0)/(2.0*gamma)
c             gamma_3 = 1.0/gamma_1
c             p_LR = (p_L/p_R)**gamma_1
c             p_star = ((cs_L + cs_R - 0.5*(gamma - 1.0)*(u_R - u_L))
c     &            /(cs_l/(p_L**gamma_1) + cs_R/(p_R**gamma_1)))**gamma_3
c             u_star = (p_LR*u_L/cs_L + u_R/cs_R + 
c     &                 2.0*(p_LR - 1.0)/(gamma - 1.0))
c     &                 /(p_LR/cs_L + 1.0/cs_R)         
c             
c           else                     !Use Two-Shock Riemann Solver (TSRS)
c             p_0 = max(0.,p_PVRS)
c             A_L = 2.0/((gamma + 1.0)*rho_L)
c             B_L = p_L*(gamma - 1.0)/(gamma + 1.0)
c             g_L = sqrt(A_L/(p_0 + B_L))
c             A_R = 2.0/((gamma + 1.0)*rho_R)
c             B_R = p_R*(gamma - 1.0)/(gamma + 1.0)
c             g_R = sqrt(A_R/(p_0 + B_R))
c             p_star = (g_L*p_L + g_R*p_R - (u_R - u_L))/(g_L + g_R)
c             u_star = 0.5*(u_R + u_L) 
c     &              + 0.5*((p_star - p_R)*g_R - (p_star - p_L)*g_L)
c           endif
c         end if
c        end if  !End of:  if(1.eq.2)then
      end if

      
      !--- Left Shock Speed ---
      if(p_star.lt.p_L)then      !Rarefaction
        q_L = 1.0
      else if(p_star.gt.1.00000*p_L)then         !Shock  1.00001*p_L: 1.00001 is needed for some gfortran bug not recognising p_star.gt.p_L !!!!!
        if(i_EoS.eq.2)then  !Euler equations
          rho_starL = rho_L*((p_star/p_L)**(1.0/gamma))
          q_L = sqrt(((p_L - p_star)/((rho_L/rho_starL) - 1.))
     &        /(gamma*(p_L)))
        else    !Tait EoS for nearly incompressible water
          d_rho_starL = rho_L - rho_star
          if(d_rho_starL.gt.0.0)then
            q_L = sqrt((rho_star/rho_L)*(p_star - p_L)/(d_rho_starL))
     &            /cs_L
          else
            q_L = 1.0
          endif
        end if
      else
        q_L = 1.0
      end if
      S_L = u_L - cs_L*q_L

      !--- Right Shock Speed ---
      if(p_star.lt.p_R)then   !Rarefaction
        q_R = 1.0
      else  if(p_star.gt.1.00000*p_R)then     !Shock   !1.00001*p_R: 1.00001 is needed for some gfortran bug not recognising p_star.gt.p_L !!!!!
        if(i_EoS.eq.2)then  !Euler equations
          rho_starR = rho_R*((p_star/p_R)**(1.0/gamma))
          q_R = sqrt(((p_R - p_star)/((rho_R/rho_starR) - 1.))
     &        /(gamma*(p_R)))
        else
          d_rho_starR = rho_R - rho_star
          if(d_rho_starR.gt.0.0)then
            q_R = sqrt((rho_star/rho_R)*(p_star - p_R)/(d_rho_starR))
     &            /cs_R
          else
            q_R = 1.0
          endif
        end if
      else
        q_R = 1.0
      end if
      S_R = u_R + cs_R*q_R
      
      !--- Contact discontinuity Shock Speed ---
      S_star = u_star
      
      
      !x_over_t = u_0   !0.0
      x_over_t = u_0_R
      if(S_L.gt.x_over_t)then
        rho_HLLC    = rho_L 
        u_HLLC      = u_L
c       v_HLLC      = v_L
        w_HLLC      = w_L
        E_HLLC      = E_L      
      elseif(S_star.ge.x_over_t)then    !elseif(S_L   .le.x_over_t.and.S_star.ge.x_over_t)then
        Q_factor_L  = rho_L*((S_L - u_L)/(S_L - S_star))
        rho_HLLC    = Q_factor_L
        rho_u_HLLC  = S_star*Q_factor_L
c       rho_v_HLLC  = v_L   *Q_factor_L
        rho_w_HLLC  = w_L   *Q_factor_L
        E_HLLC      = Q_factor_L*((E_L/rho_L) + 
     &               (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))))
        u_HLLC = rho_u_HLLC/rho_HLLC
c       v_HLLC = rho_v_HLLC/rho_HLLC
        w_HLLC = rho_w_HLLC/rho_HLLC      
      elseif(S_R   .ge.x_over_t)then    !elseif(S_star.le.x_over_t.and.S_R   .ge.x_over_t)then
        Q_factor_R  = rho_R*((S_R - u_R)/(S_R - S_star))
        rho_HLLC    = Q_factor_R
        rho_u_HLLC  = S_star*Q_factor_R
c       rho_v_HLLC  = v_R   *Q_factor_R
        rho_w_HLLC  = w_R   *Q_factor_R
        E_HLLC      = Q_factor_R*((E_R/rho_R) + 
     &               (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))))
        u_HLLC = rho_u_HLLC/rho_HLLC
c       v_HLLC = rho_v_HLLC/rho_HLLC
        w_HLLC = rho_w_HLLC/rho_HLLC
      elseif(S_R   .lt.x_over_t)then    !elseif(S_R   .le.x_over_t)then
        rho_HLLC    = rho_R 
        u_HLLC      = u_R
c       v_HLLC      = v_R
        w_HLLC      = w_R
        E_HLLC      = E_R
      else
        print*
        print*,'Problem in Riemann Solver'
        print*,'itime ',itime
        print*,'i, j ',i,j
        print*,'xp(i), zp(i) ',xp(i), zp(i)
        print*,'xp(j), zp(j) ',xp(j), zp(j)
        print*,'up(i), wp(i) ',up(i), wp(i)
        print*,'up(j), wp(j) ',up(j), wp(j)
        print*,'u_L, w_L ',u_L, w_L
        print*,'u_R, w_R ',u_R, w_R
        print*,'rho_L, p_L ',rho_L, p_L
        print*,'rho_R, p_R ',rho_R, p_R
        print*,'cs_L, cs_R ',cs_L, cs_R
        print*,'E_L, E_R ',E_L, E_R
        print*,'p_L, p_R, p_star ',p_L, p_R, p_star
        print*,'q_L, q_R ',q_L, q_R
        !print*,'rho_starL, rho_starR ',rho_starL, rho_starR
        print*,'d_rho_starL, d_rho_starR ',d_rho_starL, d_rho_starR
        print*,'rho_star, cs_star ',rho_star, cs_star
        print*,'S_L, S_star, S_R ',S_L, S_star, S_R
        print*,'p_star.gt.p_L ',p_star.gt.p_L
        print*,'p_star.gt.p_R ',p_star.gt.p_R
        print*,'(p_R*rho_L*cs_L + p_L*rho_R*cs_R - 
     &rho_L*cs_L*rho_R*cs_R*(u_R - u_L))
     &          /(rho_L*cs_L + rho_R*cs_R)'
        write(*,'(e24.16)')(p_R*rho_L*cs_L + p_L*rho_R*cs_R - 
     &                 rho_L*cs_L*rho_R*cs_R*(u_R - u_L))
     &                /(rho_L*cs_L + rho_R*cs_R)
        write(*,800)'p_L    ',p_L
        write(*,800)'rho_L  ',rho_L
        write(*,800)'cs_L   ',cs_L
        write(*,800)'u_L    ',u_L
        write(*,800)'p_R    ',p_R
        write(*,800)'rho_R  ',rho_R
        write(*,800)'cs_R   ',cs_R
        write(*,800)'u_R    ',u_R
        write(*,800)'p_star ',p_star
        stop
      end if
800   format(a8,' ',e24.16)        
      
      rho_a    =  rho_HLLC
      u_a      =  u_HLLC*cos_theta - w_HLLC*sin_theta
      w_a      =  u_HLLC*sin_theta + w_HLLC*cos_theta

      E_a      =  E_HLLC
      TE_a     =  E_a/rho_a - 0.5*(u_a*u_a + w_a*w_a)
      call equation_of_state(rho_a,TE_a,p_a,cs_a)

      return
      end
