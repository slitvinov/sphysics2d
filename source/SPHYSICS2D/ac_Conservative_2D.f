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

      subroutine ac_main
c
      include 'common.2D'
            

c
c  ...  store useful arrays
c
      !- Need to zero each object for multiobjects -
      bigUdot   = 0.0
      bigWdot   = 0.0
        
      X_Friction   = 0.0
      Z_Friction   = 0.0
        
      nb_inFriction = 0

      do i=nbfm+1,nb

c        -- Zeroing Variables for Free-Moving Objects --         

         ax(i) = 0.
         az(i) = 0.

         ar(i) = 0.
         aVol(i) = 0.
         aTE(i) = 0.

      enddo

      do i=nstart,np

c        -- Zeroing Variables --         

         ax(i) = 0.
         az(i) = 0.

         ar(i) = 0.
         aVol(i) = 0.
         aTE(i) = 0.

         ux(i) = 0.
         wx(i) = 0.

         
         sum_wab(i) = 0.
         rhop_sum(i) = 0.

         dudx_CSPHo(i) = dudx_CSPH(i)
         dudz_CSPHo(i) = dudz_CSPH(i)
         dwdx_CSPHo(i) = dwdx_CSPH(i)
         dwdz_CSPHo(i) = dwdz_CSPH(i)
         dTEdx_CSPHo(i) = dTEdx_CSPH(i)
         dTEdz_CSPHo(i) = dTEdz_CSPH(i)
         drhodx_CSPHo(i) = drhodx_CSPH(i)
         drhodz_CSPHo(i) = drhodz_CSPH(i)
         
         drhoudx_CSPHo(i) = drhoudx_CSPH(i)
         drhoudz_CSPHo(i) = drhoudz_CSPH(i)
         drhowdx_CSPHo(i) = drhowdx_CSPH(i)
         drhowdz_CSPHo(i) = drhowdz_CSPH(i)
         dTotalEdx_CSPHo(i) = dTotalEdx_CSPH(i)
         dTotalEdz_CSPHo(i) = dTotalEdz_CSPH(i)

         dudx_CSPH(i) = 0.0
         dudz_CSPH(i) = 0.0
         dwdx_CSPH(i) = 0.0
         dwdz_CSPH(i) = 0.0
         dTEdx_CSPH(i) = 0.0
         dTEdz_CSPH(i) = 0.0
         drhodx_CSPH(i) = 0.0
         drhodz_CSPH(i) = 0.0
         
         drhoudx_CSPH(i) = 0.0
         drhoudz_CSPH(i) = 0.0
         drhowdx_CSPH(i) = 0.0
         drhowdz_CSPH(i) = 0.0
         dTotalEdx_CSPH(i) = 0.0
         dTotalEdz_CSPH(i) = 0.0

c       ---  Clear vorticity arrays ---

         if(ipoute.eq.1.and.i_vort.eq.1)then         
           curl3y_a(i)= 0.
           curl3y_b(i)= 0.
         end if

      enddo

       ly2 = 0   !Default Y-Periodic value
       !kind=1 !kind1+kind2-1 
       do kind_p1=1,2
       ini_kind_p2=mod(kind_p1,2)+1
       
       
       do lz=1,ncz
       do lx=1,ncx
         j1 = lx + (lz-1)*ncx
         
c         print*,'kind_p1, lx, lz, j1 ',kind_p1,lx, lz, j1
         
         if(nc(j1,kind_p1).gt.0) then
c                            ! if the cell is not empty, then
c                            ! loop over it and over neighboring
c                            ! cells


c          -- Cells in the same XY sheet -
           lx2 = lx+1
           if(lx2.le.ncx)then
             call celij(j1,j1+1,kind_p1,ini_kind_p2,lx2)     !East
           endif      
      
           lz2=lz+1
           if(lz2.le.ncz)then
       
             !- Same row -
             call celij(j1,j1+ncx,kind_p1,ini_kind_p2,lx2)   !North
       
             lx2=lx+1
             if(lx2.gt.2) call celij(j1,j1+ncx-1,
     &                                 kind_p1,ini_kind_p2,lx2)   !North-West
       
             lx2=lx+1
             if(lx2.le.ncx) call celij(j1,j1+ncx+1,
     &                                 kind_p1,ini_kind_p2,lx2)   !North-East                 
       
           endif    !End of:   if(lz2.le.ncz)then
         endif    !End of:  if(nc(j1,kind_p1).gt.0) then
       
       enddo
       enddo
       
       do j1=1,nct
           if(nc(j1,kind_p1).gt.0)
     +		call self(j1,kind_p1,ini_kind_p2)
       enddo
             
c      -- Periodic Boundary Calls in X-Direction --              
       if(i_periodicOBs(1).eq.1)then
         !- Special Treatment for Right Column Cells (lx=ncx) -
         lx = ncx
         lx2 = lx+1
         do lz = 1,ncz-1
           j1 = lx + (lz-1)*ncx
           call celij(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
           call celij(j1,j1+ncall13,kind_p1,ini_kind_p2,lx2) !North-East X-Periodic
         end do    !End of:  do lz = 1,ncz-1       
         !- Special Treatment for Corner Cell (lx=ncx, lz=ncz) -
         !- Note this is for X-Periodicity Only!               -
         !Norths & North-Easts of lx = 1 
         lz = ncz
         j1 = lx + (lz-1)*ncx
         call celij(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
         !- Special Treatment for Left Column Cells (lx=1) -
         lx = 1
         lx2 = lx+1
         do lz = 1,ncz-1
           j1 = lx + (lz-1)*ncx
           call celij(j1,j1+ncall12,kind_p1,ini_kind_p2,lx2) !North-West X-Periodic
         end do    !End of:  do lz = 1,ncz-1       
       endif   !End of:  if(i_periodicOBs(1).eq.1)then
       !-----  End of Y-Periodicity --------------------------   
             
       enddo   !End of:  do kind_p1=1,2
             
       
       do i=nstart,np
         
c         udot(i) = ax(i)
c         wdot(i) = az(i)
c         rdot(i) = ar(i)
       
         if(sum_wab(i).lt.0.1) then
           one_over_sum_wab=1.
         else
           one_over_sum_wab=1./sum_wab(i)
         endif
       
         xcor(i) = eps*ux(i)   !*one_over_sum_wab	! NORMALIZING KERNELS
         zcor(i) = eps*wx(i)   !*one_over_sum_wab
       
c         TEdot(i)=aTE(i)
       
       enddo
              
       return
       end

