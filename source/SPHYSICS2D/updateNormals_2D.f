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

      subroutine updateNormals_2D(i_start, i_stop)

      include 'common.2D'
      
            
    
c      i_start = 1         !BedPs only
c      i_stop  = nbf       !Bed particles only
c
c      i_start = nbf + 1   !Moving Solid Object particles only
c      i_stop  = nb        !Moving Solid Object particles only

      !-- Serial Version --
      do i = i_start,i_stop  
        xp_i = xp(i)
        zp_i = zp(i)
        
        iminus1_neighbour = iBP_Pointer_Info(i,3)
        iplus1_neighbour  = iBP_Pointer_Info(i,4)

        !xp_BP_iminus1 = BP_xz_Data(iminus1_neighbour,1)
        !zp_BP_iminus1 = BP_xz_Data(iminus1_neighbour,2)
        !xp_BP_iplus1 = BP_xz_Data(iplus1_neighbour,1)
        !zp_BP_iplus1 = BP_xz_Data(iplus1_neighbour,2)
        xp_BP_iminus1 = xp(iminus1_neighbour)
        zp_BP_iminus1 = zp(iminus1_neighbour)
        xp_BP_iplus1 = xp(iplus1_neighbour)
        zp_BP_iplus1 = zp(iplus1_neighbour)
        
        if(i.eq.iminus1_neighbour)then
          xp_previous = xp_i - (xp_BP_iplus1-xp_i)
          zp_previous = zp_i - (zp_BP_iplus1-zp_i)
          xp_next = xp_BP_iplus1
          zp_next = zp_BP_iplus1
        else if(i.eq.iplus1_neighbour)then
          xp_previous = xp_BP_iminus1
          zp_previous = zp_BP_iminus1
          xp_next = xp_i + (xp_i-xp_BP_iminus1)
          zp_next = zp_i + (zp_i-zp_BP_iminus1)
        else
          xp_previous = xp_BP_iminus1
          zp_previous = zp_BP_iminus1
          xp_next = xp_BP_iplus1
          zp_next = zp_BP_iplus1
        end if
        
        dxn_forward=(xp_next-xp_i)
        dzn_forward=(zp_next-zp_i)
        dxn_back=(xp_i-xp_previous)
        dzn_back=(zp_i-zp_previous)
        

        !- X-Periodic Correction -
        if(i_periodicOBs(1).eq.1.and.i.gt.nbfm)then
          
          !- Tangent t Check Fwd -
          drx_check = dxn_forward
          xp_i_check = xp_next
          xp_j_check = xp_i
          if(drx_check.gt.four_h)then
            dxp_dble1 = xmax_container_double-dble(xp_i_check)
            dxp_dble2 = dble(xp_j_check)-xmin_container_double
            drx_check = -real(dxp_dble1+dxp_dble2)
          end if
          if(-drx_check.gt.four_h)then
            dxp_dble1 = dble(xp_i_check)-xmin_container_double
            dxp_dble2 = xmax_container_double-dble(xp_j_check)
            drx_check =  real(dxp_dble1+dxp_dble2)
          end if
          dxn_forward = drx_check
          
          !- Tangent s Check Backwd -
          drx_check = dxn_back
          xp_i_check = xp_i
          xp_j_check = xp_previous
          if(drx_check.gt.four_h)then
            dxp_dble1 = xmax_container_double-dble(xp_i_check)
            dxp_dble2 = dble(xp_j_check)-xmin_container_double
            drx_check = -real(dxp_dble1+dxp_dble2)
          end if
          if(-drx_check.gt.four_h)then
            dxp_dble1 = dble(xp_i_check)-xmin_container_double
            dxp_dble2 = xmax_container_double-dble(xp_j_check)
            drx_check =  real(dxp_dble1+dxp_dble2)
          end if
          dxn_back = drx_check
        
        endif
        
        rbn_back = sqrt(dxn_back*dxn_back + dzn_back*dzn_back)
        rbn_forward = sqrt(dxn_forward*dxn_forward +
     &                     dzn_forward*dzn_forward)
        dzn = dzn_forward + dzn_back
        dxn = dxn_forward + dxn_back
        rb=dxn*dxn+dzn*dzn
        xnb(i)=-dzn/sqrt(rb)
        znb(i)=dxn/sqrt(rb)
        xtb(i) =  znb(i)
        ztb(i) = -xnb(i)
        if(abs(rbn_back).gt.0.0)then
          deltaptb(i,1) = rbn_back
        else
          deltaptb(i,1) = h*1.0e-3
        end if          
        if(abs(rbn_forward).gt.0.0)then
          deltaptb(i,2) = rbn_forward
        else
          deltaptb(i,2) = h*1.0e-3
        end if 
         
         
      enddo

      end
