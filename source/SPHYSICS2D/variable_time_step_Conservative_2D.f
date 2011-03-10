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

      subroutine variable_time_step
      include 'common.2D'
      

C	Calculate time step 
c	Conditions given in 
c      Monaghan (1989), JCP 82, 1-15
c      Monaghan & Koss (1999), JWPCOE 125, 145-154

      Fa_max=-1000.0
      cs_max=-1000.0
      i_fa = -1000

      if(iBC.eq.1)then
        nstart_VTS = nstart
      else
        nstart_VTS = nbp1
      endif

      do i=nstart_VTS,np
        Fa=sqrt((ax(i) - Volrho(i)*grx)**2
     &        + (az(i) - Volrho(i)*grz)**2)

        if (Fa.ge.Fa_max) then
           i_fa=i
           Fa_max=Fa
        endif
        if (cs(i).ge.cs_max) then
           cs_max=cs(i)
        endif
      enddo 

      if(num_FB.gt.0)then
      !print*,'itime, dt ',itime,dt
      !print*,'i_fa, Fa_max ',i_fa, Fa_max
      do i=nbfm+1,nb
        i_num_FB = i_FB_Pointer_Info(i)
        Fa=sqrt((ax(i) )**2
     &        + (az(i) )**2)

        if (Fa.ge.Fa_max) then
           i_fa=i
           Fa_max=Fa
        endif
        if(cs(i).ge.cs_max)then
          cs_max=cs(i)
        endif
        vel_magnitude = sqrt(up(i)*up(i) + wp(i)*wp(i))
        if (vel_magnitude.ge.cs_max) then
           cs_max=vel_magnitude
        endif
      enddo
      if(.NOT.(abs(Fa_max).gt.0.0))then
        i_fa = 0
      endif
      !print*,'i_fa, Fa_max ',i_fa, Fa_max
      endif    !End of:   if(num_FB.gt.0)then

c      Fa_sqrt=((ax(i_fa)+grx)**2+
c     +         (az(i_fa)+grz)**2)**0.25
c      dt_1=sqrt(h)/Fa_sqrt
      if(i_fa.gt.0)then
        Fa_sqrt=((ax(i_fa)**2
     &           +az(i_fa)**2)**0.25)
     &          /sqrt(Volrho(i_fa))
        dt_1=sqrt(h)/Fa_sqrt
      else
        Fa_sqrt=((grx)**2+(grz)**2)**0.25
        dt_1=sqrt(h)/Fa_sqrt
      endif

      !dt_2=h/(cs_max+h*visc_dt)
      if(cs_max.gt.0.0)then      
        dt_2=h/(cs_max+h*visc_dt)       
      else       
        dt_2=dt_1         
      endif

      dt_new=CFL_number*min(dt_1,dt_2)
       

      return

      end
