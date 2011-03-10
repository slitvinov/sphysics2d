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
       
       program SPHYSICS
       
C
c			Axis
c			^
c			|      / Y
c			| Z   /
c			|    /
c			|   /
c			|  /
c			| /
c			|/
c			--------------> X
c
       
       include 'common.2D'
      
       
       character supp*4,name*40,detsupp*4
              
       REAL time_begin, time_end              
       CALL CPU_TIME (time_begin)
       
       open(80,file='sph.out')
       
       write(80,*) ' '
       write(80,*) '<SPHYSICS>  Copyright (C) <2009>'
       write(80,*) 
     & '<Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, '
       write(80,*) 
     & 'Dr Benedict Rogers, Dr Alejandro Crespo, '
       write(80,*) 
     & 'Dr Muthukumar Narayanaswamy, Dr Shan Zou, & Dr Andrea Panizzo >'
       write(80,*) 'This program comes with ABSOLUTELY NO WARRANTY;    '
       write(80,*) 'This is free software, and you are welcome to      '
       write(80,*) 'redistribute it under conditions stated in         '
       write(80,*) 'the GPL License;                                   '

       write(80,*) ' '
       write(80,*)' ---     SPHYSICS_2D.F       ---'
       write(80,*)' ---   Distributed under     ---'
       write(80,*)' ---    the GPL License      ---'
       write(80,*)
c      -- screen printout ---              
       print*
       write(*,*) '<SPHYSICS>  Copyright (C) <2009>'
       write(*,*) 
     & '<Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, '
       write(*,*) 
     & 'Dr Benedict Rogers, Dr Alejandro Crespo, '
       write(*,*) 
     &'Dr Muthukumar Narayanaswamy, Dr Shan Zou, & Dr Andrea Panizzo >'
       write(*,*) 'This program comes with ABSOLUTELY NO WARRANTY;    '
       write(*,*) 'This is free software, and you are welcome to      '
       write(*,*) 'redistribute it under conditions stated in         '
       write(*,*) 'the GPL License;                                   '

       write(*,*) ' '
       write(*,*) ' ---     SPHYSICS_2D.F       ---'
       write(*,*) ' ---   Distributed under     ---'
       write(*,*) ' ---    the GPL License      ---'
       write(*,*) ' '
       
c
c  ...  file specifications
c
       
       
       open(19,file='DT')       
       open(50,file='ENERGY')
       open(54,file='MovingBodyHistory.out')

       
       pi=4.*atan(1.)
       g=9.81
       
c
c  ...  initialization
c      
       call getdata
       call energy(g)
       
       call ini_divide(1)
       call divide(1,nbf,1)
       call keep_list
       
       if(i_restartRun.eq.1.or.i_restartRun.eq.3) then
         open(19,file='DT',status='old',access='append')
       else
         open(19,file='DT') 
         write(19,*) 'time dt1 dt2 dt_new'
       endif
       close(19)
       
c
c  .............................. MAIN LOOP ...................................
c
       
       tdetail=0
       ngrabdet=0
       grab_P=0
       grab_E=0
       grx=0.0
       gry=0.0
       grz=-g
                           
       do while (time.lt.tmax)
          visc_dt=0.
       
	    if((time-trec_ini)-out*ngrab.ge.out) then 
			ipoute =1
		else
			ipoute =0
		endif
	   
	              	         
         !-- Call Main Subroutines --  
         call step
              
         itime=itime+1
         time=time+dt 
	      
         !-- File Output Routines --
	   if(ipoute.eq.1) then  
 
           ngrab=ngrab+1
           write(supp,'(i4.4)') ngrab
           name='PART_'//supp
           write(80,*) name
           write(*,*) name
           open(23,file=name)

           if(i_vort.eq.1) then
               name='VORT_'//supp
               open(24,file=name)
           endif   
		     
           call poute(23)
           
           close(23)
           close(24)
       
           open(44,file='RESTART')
           write(44,100)itime,time,ngrab,dt  
           close(44)
           
           grab_P=0.
       
         else

           grab_P=grab_P+dt       
         endif


100   format(i8,e16.8,i8,e16.8) 
       
         
         if(i_densityFilter.gt.0.and.
     &        itime/ndt_FilterPerform*ndt_FilterPerform.eq.itime)then 
            call densityFilter
         endif
         
         if(grab_E.gt.out.or.time.eq.trec_ini) then
           call energy(g)
           grab_E=0.
         else
           grab_E=grab_E+dt
         endif   

                                   
c        -- Detailed recording --    
          tdetail=tdetail+dt
       
          if(time.ge.t_sta_det.and.
     +       time.le.t_end_det.and.
     +       tdetail.ge.dtrec_det)then
            tdetail=0
	      ngrabdet=ngrabdet+1
	      write(detsupp,'(i4.4)') ngrabdet
            name='DETPART_'//detsupp
            write(*,*) name
            open(53,file=name)
            call poute(53)
            close(53)
          endif

                   
       enddo
       
       open(33,file='EPART')                                             
       call poute(33)   !  Final record   
       close(33)                                             
       
       
C	Closing files

c       close(19)
       close(50)       
       close(53) 
       close(54)

	   
       write(80,*) 'End '
       write(*,*) 'End '      
       
       CALL CPU_TIME ( time_end )
       write(80,*)'time_begin',time_begin,'seconds'
       write(80,*)'time_end',time_end,'seconds'
       write(80,*)'Time of operation was',time_end-time_begin,'seconds'
       write(*,*)'time_begin',time_begin,'seconds'
       write(*,*)'time_end',time_end,'seconds'
       write(*,*)'Time of operation was',time_end-time_begin,'seconds'
       
       close(80)      
       
       end       !end main program
       
