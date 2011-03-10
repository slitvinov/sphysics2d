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

      subroutine ac_MLS
c
      include 'common.2D'

c
c  ...  store useful arrays
c
            
      ly2 = 0   !Default Y-Periodic value

      do kind_p1=1,2
       ini_kind_p2=mod(kind_p1,2)+1
       
       
       do lz=1,ncz
       do lx=1,ncx
         j1 = lx + (lz-1)*ncx
         
         if(nc(j1,kind_p1).gt.0) then
c                            ! if the cell is not empty, then
c                            ! loop over it and over neighboring
c                            ! cells

c          -- Cells in the same XY sheet -
           lx2 = lx+1
           if(lx2.le.ncx)then
             call pre_celij_MLS(j1,j1+1,kind_p1,ini_kind_p2,lx2)     !East
           endif      
           lz2=lz+1
           if(lz2.le.ncz)then       
             !- Same row -
             call pre_celij_MLS(j1,j1+ncx,kind_p1,ini_kind_p2,lx2)   !North       
             lx2=lx+1
             if(lx2.gt.2) call pre_celij_MLS(j1,j1+ncx-1,
     &                                 kind_p1,ini_kind_p2,lx2)   !North-West       
             lx2=lx+1
             if(lx2.le.ncx) call pre_celij_MLS(j1,j1+ncx+1,
     &                                 kind_p1,ini_kind_p2,lx2)   !North-East   
            
                            
           endif    !End of:   if(lz2.le.ncz)then
           
         endif    !End of:  if(nc(j1,kind_p1).gt.0) then       
       
       enddo
       enddo
       
       do j1=1,nct
           if(nc(j1,kind_p1).gt.0)
     +        call pre_self_MLS(j1,kind_p1,ini_kind_p2)
       enddo
      
c      -- Periodic Boundary Calls in X-Direction --              
       if(i_periodicOBs(1).eq.1)then
         !- Special Treatment for Right Column Cells (lx=ncx) -
         lx = ncx
         lx2 = lx+1
         do lz = 1,ncz-1
           j1 = lx + (lz-1)*ncx
           call pre_celij_MLS(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
           call pre_celij_MLS(j1,j1+ncall13,kind_p1,ini_kind_p2,lx2) !North-East X-Periodic
         end do    !End of:  do lz = 1,ncz-1       
         !- Special Treatment for Corner Cell (lx=ncx, lz=ncz) -
         !- Note this is for X-Periodicity Only!               -
         !Norths & North-Easts of lx = 1 
         lz = ncz
         j1 = lx + (lz-1)*ncx
         call pre_celij_MLS(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
         !- Special Treatment for Left Column Cells (lx=1) -
         lx = 1
         lx2 = lx+1
         do lz = 1,ncz-1
           j1 = lx + (lz-1)*ncx
           call pre_celij_MLS(j1,j1+ncall12,kind_p1,ini_kind_p2,lx2) !North-West X-Periodic
         end do    !End of:  do lz = 1,ncz-1       
       endif   !End of:  if(i_periodicOBs(1).eq.1)then
       !-----  End of Y-Periodicity --------------------------   
             
      enddo   !End of:  do kind_p1=1,2  
       
      end subroutine
