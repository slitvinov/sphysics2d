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
       
      subroutine divide(n_start,n_end,kind_p)
      include 'common.2D'
      

c
       do k=n_start,n_end
         if (iflag(k).ne.0) then
              
           dx = xp(k) - xmin
           dz = zp(k) - zmin
       
           icell = int( dx * one_over_2h ) + 1
           kcell = int( dz * one_over_2h ) + 1
       
           if(i_periodicOBs(1).eq.1.and.icell.gt.ncx)icell=ncx   !Periodic Open Boundary correction X-Direction   
           !if(i_periodicOBs(2).eq.1.and.jcell.gt.ncy)jcell=ncy    !Periodic Open Boundary correction Y-Direction   
           !if(i_periodicOBs(3).eq.1.and.kcell.gt.ncz)kcell=ncz   !Periodic Open Boundary correction Z-Direction   

           ! ii is the linear cell position in the matrix of 2h cells
           ii    = icell + (kcell - 1)*ncx 

          if(ii.lt.1)then
           write(80,*)' '
           write(80,*)'ERROR in divide_2D.f'
           write(80,*)'ii.lt.1'
           write(80,*)'itime',itime
           write(80,*)'kind_p       ',kind_p
           write(80,*)'nplink_max',nplink_max
           write(80,*)'n_start,n_end ',n_start,n_end
           write(80,*)'k',k
           write(80,*)'xp(k), zp(k)',xp(k), zp(k)
           write(80,*)'icell, kcell ',icell, kcell
           write(80,*)'Box ii ',ii
           write(80,*)'nc(ii,kind_p)',nc(ii,kind_p)
c          -- screen printout ---              
           write(*,*)' '
           write(*,*)'ERROR in divide_2D.f'
           write(*,*)'ii.lt.1'
           write(*,*)'itime',itime
           write(*,*)'kind_p       ',kind_p
           write(*,*)'nplink_max',nplink_max
           write(*,*)'n_start,n_end ',n_start,n_end
           write(*,*)'k',k
           write(*,*)'xp(k), zp(k)',xp(k), zp(k)
           write(*,*)'icell, kcell ',icell, kcell
           write(*,*)'Box ii ',ii
           write(*,*)'nc(ii,kind_p)',nc(ii,kind_p)
           stop
          elseif(ii.gt.nct)then
           write(80,*)' '
           write(80,*)'ERROR in divide_2D.f'
           write(80,*)'ii.lt.1'
           write(80,*)'itime',itime
           write(80,*)'kind_p       ',kind_p
           write(80,*)'nplink_max',nplink_max
           write(80,*)'n_start,n_end ',n_start,n_end
           write(80,*)'k',k
           write(80,*)'xp(k), zp(k)',xp(k), zp(k)
           write(80,*)'icell, kcell ',icell, kcell
           write(80,*)'Box ii ',ii
           write(80,*)'nc(ii,kind_p)',nc(ii,kind_p)
c          -- screen printout ---              
           write(*,*)' '
           write(*,*)'ERROR in divide_2D.f'
           write(*,*)'ii.lt.1'
           write(*,*)'itime',itime
           write(*,*)'kind_p       ',kind_p
           write(*,*)'nplink_max',nplink_max
           write(*,*)'n_start,n_end ',n_start,n_end
           write(*,*)'k',k
           write(*,*)'xp(k), zp(k)',xp(k), zp(k)
           write(*,*)'icell, kcell ',icell, kcell
           write(*,*)'Box ii ',ii
           write(*,*)'nc(ii,kind_p)',nc(ii,kind_p)
           stop
          end if
                                                        
           ! nc is the number of particles in cell ii
           nc(ii,kind_p) = nc(ii,kind_p)+1       
                                      
          if(nc(ii,kind_p).gt.nplink_max)then
           write(80,*) ' '
           write(80,*)'ERROR in divide_2D.f'
           write(80,*)'nc(ii,kind_p) >= nplink_max'
           write(80,*)'itime',itime
           write(80,*)'nc(ii,kind_p)',nc(ii,kind_p)
           write(80,*)'kind_p       ',kind_p
           write(80,*)'nplink_max',nplink_max
           write(80,*)'k',k
           write(80,*)'xp(k), zp(k)',xp(k), zp(k)
           write(80,*)'icell, kcell ',icell, kcell
           write(80,*)'Box ii ',ii
c          -- screen printout ---              
           write(*,*) ' '
           write(*,*)'ERROR in divide_2D.f'
           write(*,*)'nc(ii,kind_p) >= nplink_max'
           write(*,*)'itime',itime
           write(*,*)'nc(ii,kind_p)',nc(ii,kind_p)
           write(*,*)'kind_p       ',kind_p
           write(*,*)'nplink_max',nplink_max
           write(*,*)'k',k
           write(*,*)'xp(k), zp(k)',xp(k), zp(k)
           write(*,*)'icell, kcell ',icell, kcell
           write(*,*)'Box ii ',ii
           stop
          end if
          
          ibox(ii,kind_p,nc(ii,kind_p))=k  !Tells us that particle with array location k (i.e. xp(k) )  
                                         !is in box ii which, so far, contains "nc(ii,mkind)" particles
           
         endif
         
c         if(k.eq.39650)then
c         if(k.eq.1093.or.k.eq.1103.or.
c     &      k.eq.39650.or.k.eq.39700.or.k.eq.57574)then
c           print*
c           print*,'in DIVIDE'
c           print*,'k',k
c           print*,'xp(k), zp(k)',xp(k),zp(k)
c           print*,'icell,kcell',icell,kcell
c           print*,'box ii',ii
c           print*,'nc(ii,kind_p) ',nc(ii,kind_p)
c           print*,'ncx , ncz  ',ncx, ncz 
c           print*,'kind_p ',kind_p
c           print*,'n_start, n_end ',n_start, n_end
c           print*,'one_over_2h ',one_over_2h
c           print*,'dx, dz ',dx,dz
c           print*,'dx * one_over_2h  ',dx * one_over_2h
c           print*,'int(dx * one_over_2h)  ',int(dx * one_over_2h)
c           print*
c           read(*,*)
c         end if

       enddo 
       
c      print*,'DIVIDE itime, tmax,k,n_start,n_end ',
c     &               itime,tmax,k,n_start,n_end  

c      itime_check = 0
c      if(itime.ge.itime_check)then   !.and.mkind.eq.1)then
c        print*
c        print*,'At end of Divide, kind_p = ',kind_p
c        print*,'itime ',itime
c        print*,'n_start, n_end ',n_start,n_end
c        !print*,'xmin_local, ymin_local',xmin_local, ymin_local
c        print*,'nct_local ',nct_local
c        ii = 289
c        print*,'Box ii ',ii
c        print*,'nc(ii,mkind=1) ',nc(ii,1)
c        print*,'nc(ii,mkind=2) ',nc(ii,2)
c        print*,'nc_ii_max',nc_ii_max
c        !print*
c          print*,'For each kind'
c          do mmkind = 1,2
c            print*,'mmkind ',mmkind
c            do iii = 1,nc(ii,mmkind)
c              i = ibox(ii,mmkind,iii)
c              print*,'iii, i,xp, zp ',iii,i,xp(i),zp(i)
c            end do
c            !is1 = ist(ii,mmkind)
c            !ie1 = ien(ii,mmkind)
c            !print*,'is1, ie1 ',is1,ie1
c            !if(nc(ii,mmkind).gt.0)then
c            !  do iii = is1,ie1
c            !    i = ip(iii,mmkind)
c            !    print*,'iii, i,xp, zp ',iii,i,xp(i),zp(i)
c            !  end do
c            !endif
c          end do
c        print*
c       endif

c       print*,'itime, nc(11204,1+2) ',itime, nc(11204,1), nc(11204,2)
       
       return
       end
