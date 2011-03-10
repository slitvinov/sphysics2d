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

      subroutine correct
c
      include 'common.2D'
      
c
c     ...  account for body forces
c
      do i=nbp1,np

        udot(i) = udot(i) + grx*iflag(i)
        wdot(i) = wdot(i) + grz*iflag(i)

      enddo

c
c      ...  account for XSPH
c

      do i=nbp1,np

        xdot(i) = up(i) + xcor(i)
        zdot(i) = wp(i) + zcor(i)
              
      enddo


      do i=nbp1,np

c     --- SPS terms ---

      !Smag = (0.12*dpx)**2
		           
	      Prr=2.0*((dudx_CSPH(i)**2) +(dwdz_CSPH(i)**2))
     &                   +((dudz_CSPH(i)+dwdx_CSPH(i))**2)     
    
	      visc_SPS= Smag*sqrt(Prr)   
	      div_u=dudx_CSPH(i)+dwdz_CSPH(i)
	      sps_k=twoThirds*visc_SPS*div_u 
	      sps_Blin_etal = Blin_const*Prr   !2.0*0.0066*dp*Prr     !Blin et al.

	      !sps_k = twoThirds*0.18*dp*Prr    !Gotoh et al. 2001

	      tau_xx=2.0*visc_SPS*dudx_CSPH(i)-sps_k-sps_Blin_etal
	      tau_xz=visc_SPS*(dudz_CSPH(i)+dwdx_CSPH(i))
	      tau_zz=2.0*visc_SPS*dwdz_CSPH(i)-sps_k-sps_Blin_etal

	      uno_rho2=1.0/rhop(i)         !uno_rho2 = 1.0/(rhop(i)*rhop(i))

	      tau_xx_uno_rho2(i)=tau_xx*uno_rho2  
	      tau_xz_uno_rho2(i)=tau_xz*uno_rho2
	      tau_zz_uno_rho2(i)=tau_zz*uno_rho2

	      if(visc_SPS.gt.visc_SPS_max)then
	        visc_SPS_max = visc_SPS
	        i_visc_max = i
	      end if

	enddo

      return
      end
