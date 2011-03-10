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

      subroutine poute(nf)

c
      include 'common.2D'
      
      character Normalsname*40,MovingObjectsFile*40


c     --- Write out to PART-type files ---


      do i=1,np
         write(nf,100) xp(i),zp(i),up(i),wp(i),rhop(i),p(i),pm(i)

      enddo


100   format(7e16.8)


c     --- Write out to DT file ---

      open(19,file='DT',status='old',POSITION='append')
        write(19,*) time, dt_1, dt_2, dt_new
      close(19)


c     --- Write out to VORTICITY file ---

      if(i_vort.eq.1) then
         do i=1,nb
            vorty_temp=0.0  !Vorticity in an y=const plane
            write(nf+1,105) vorty_temp
         enddo

         do i=nbp1,np
           vorty_temp= curl3y_a(i)-curl3y_b(i)  !Vorticity in an y=const plane
           write(nf+1,105) vorty_temp
         enddo

      endif

105   format(1e16.8)


        MovingObjectsFile='Floating_Bodies.RESTART'
        open(nf+2,file=MovingObjectsfile,STATUS='replace')
        write(nf+2,*)iopt_FloatingBodies
        write(nf+2,*)nbfm
        do i_num_FB=1,num_FB
          write(nf+2,*)i_num_FB
          write(nf+2,*)bigMass(i_num_FB)
          write(nf+2,*)bigInertiaYY(i_num_FB)
          write(nf+2,*)XcylinderDimension(i_num_FB),
     &                 ZcylinderDimension(i_num_FB)
          write(nf+2,*)cylinderDensity(i_num_FB)
          write(nf+2,*)FB_SpecificWeight(i_num_FB)
          write(nf+2,*)friction_coeff(i_num_FB)
          write(nf+2,*)Box_XC(i_num_FB),Box_ZC(i_num_FB)
          write(nf+2,*)bigU(i_num_FB),bigW(i_num_FB),
     &                 bigOmegaY(i_num_FB)
          write(nf+2,*)nb_FB(i_num_FB)
        enddo
        close(nf+2)

c     --- Write out to sph.out file ---
 
      write(80,*) 'In poute_2D, itime = ',itime
      write(80,*) 'time, dt  ',time,dt,'  ok'
      write(*,*) 'In poute_2D, itime = ',itime
      write(*,*) 'time, dt  ',time,dt,'  ok'


      return 
      end

