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

      subroutine densityFilter

      include 'common.2D'


      itime_check = 10e8   !0    !10e8   !0    !
      i_ParticleCheck = 23197  !41805   !231   !nb_local + 1   !3299   !
      if(itime.ge.itime_check)then
        rhop_old = rhop(i_ParticleCheck)
      endif

      nstart_Shepard = nbp1   !nstart

       call ac_Shepard

      if(itime.ge.itime_check)then
        i =  i_ParticleCheck
        sum_wab_check = sum_wab(i)
        rhop_sum_check = rhop_sum(i)
      endif

       do i=nstart_Shepard,np
         sum_wab(i) = sum_wab(i)+ adh*pm(i)/rhop(i)  !self contribution
         rhop_sum(i) = rhop_sum(i)+ adh*pm(i)  !self contribution
         rhop(i) = rhop_sum(i)/sum_wab(i)
         pVol(i) = pm(i)/rhop(i)
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))
       enddo

       if(itime.ge.itime_check)then
        print*
        print*,'In shepard.f,itime = ',itime
        i =  i_ParticleCheck
        print*,'i, xp, zp ',i,xp(i),zp(i)
        print*,'rhop_old ',rhop_old
        print*,'rhop(i)  ',rhop(i)
        print*,'p(i)    ',p(i)
        print*,'pm(i)   ',pm(i)
        print*,'sum_wab_check  ',sum_wab_check 
        print*,'rhop_sum_check ',rhop_sum_check 
        print*,'adh    ',adh
        print*,'adh*pm(i)   ',adh*pm(i)
        print*,'rhop_sum_check + adh*pm(i) ',rhop_sum_check + adh*pm(i)
        print*,'sum_wab(i)  ',sum_wab(i)
        print*,'rhop_sum(i) ',rhop_sum(i)
        print*,'pVol(i) ',pVol(i)
        print*,'nb_local ',nb
        print*,'np_local ',np
        print*
        !read(*,*)
        if(itime.ge.1000)stop
       end if

       return
       end
