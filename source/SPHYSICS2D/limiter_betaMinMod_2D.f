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

      subroutine limiter(i,j,delphi_i,delphi_j,d_phi_ij,
     &                   del_phi_lim_i,del_phi_lim_j)
c
      include 'common.2D'

        !- MUSCL Limiters - Toro (2001) p. 208
        diff_factor = 0.5
        if(delphi_j.gt.0.0)then
          del_phi_lim =max(0.,min(beta_lim*delphi_i,delphi_j),
     &                        min(delphi_i,beta_lim*delphi_j))
          if(del_phi_lim.gt.d_phi_ij.and.abs(del_phi_lim).gt.0.)then
                if(delphi_i.gt.delphi_j)then
                  diff_factor_i = diff_factor
                  diff_factor_j = 0.0
                else
                  diff_factor_i = 0.0
                  diff_factor_j = diff_factor
                endif
                del_phi_lim_i = d_phi_ij*diff_factor_i
                del_phi_lim_j = d_phi_ij*diff_factor_j
          else
                del_phi_lim_i = del_phi_lim
                del_phi_lim_j = del_phi_lim
          endif
        else
          del_phi_lim =min(0.,max(beta_lim*delphi_i,delphi_j),
     &                        max(delphi_i,beta_lim*delphi_j))
          if(del_phi_lim.lt.d_phi_ij.and.abs(del_phi_lim).gt.0.)then
                if(delphi_i.lt.delphi_j)then
                  diff_factor_i = diff_factor
                  diff_factor_j = 0.0
                else
                  diff_factor_i = 0.0
                  diff_factor_j = diff_factor
                endif
                del_phi_lim_i = d_phi_ij*diff_factor_i
                del_phi_lim_j = d_phi_ij*diff_factor_j
          else
                del_phi_lim_i = del_phi_lim
                del_phi_lim_j = del_phi_lim
          endif
        end if

      end subroutine