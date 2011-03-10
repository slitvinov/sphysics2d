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


c       initial data for SPH models                                                                     72

      program SPHYSICSgen_2D
      include 'common.gen2D'
           
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


      print*
      print*,' ---    SPHYSICSgen_2D.F     ---'
      print*,' ---   GENERATING GEOMETRY   ---'
      print*,' ---   Distributed under     ---'
      print*,' ---    the GPL License      ---'
      print*

      write(*,*) 'Choose Starting options:   Start new RUN = 0 ' 
      write(*,*) '                       : Restart old RUN = 1 '
      write(*,*) '     with CheckPointing:   Start new RUN = 2 '
      write(*,*) '                       : Restart old RUN = 3 '
      read(*,*) i_restartRun
      write(*,*) i_restartRun

C     KERNEL
      write(*,*)'Choose KERNEL'

      write(*,*)'Gaussian        = 1'
      write(*,*)'Quadratic       = 2'
      write(*,*)'Cubic- Spline   = 3'

      write(*,*)'Quintic Wendland= 5'
      read(*,*) i_kernel
      write(*,*) i_kernel

C     ALGORITHM
      write(*,*)'Choose ALGORITHM'
      write(*,*)'Predictor- Corrector = 1'
      write(*,*)'Verlet               = 2'
      write(*,*)'Symplectic           = 3'
      write(*,*)'Beeman               = 4'
      read(*,*) i_algorithm
      write(*,*) i_algorithm

      if(i_algorithm.eq.2)then
        write(*,*)' ndt_VerletPerform ? '
        read(*,*) ndt_VerletPerform
        write(*,*) ndt_VerletPerform
      else
	  ndt_VerletPerform=0
      endif
	  

c     Density filter

      write(*,*) ' Use Density filter (0=no, 1=Shepard, 2=MLS) ??'
      read(*,*) i_densityFilter
      write(*,*) i_densityFilter

      if(i_densityFilter.gt.0)then
        write(*,*)' ndt_FilterPerform ? '
        read(*,*) ndt_FilterPerform
        write(*,*) ndt_FilterPerform
      else
        ndt_FilterPerform=0
      endif


c     Kernel & Gradient corrections
      write(*,*) ' Use kernel & gradient correction ??'
      write(*,*) ' 0=no, 1= Kernel correction, 
     + 2= Gradient kernel correction'
      read(*,*)  i_kernelcorrection
      write(*,*) i_kernelcorrection

c     Constants
c
      g=9.81
      pi=4.0*atan(1.0)


c     Initialise periodicity
c
      i_periodicOBs(1) = 0
      i_periodicOBs(2) = 0
      i_periodicOBs(3) = 0



c	Viscosity terms
      write(*,*) 'Kind of viscosity  ??'
      write(*,*) ' (1) Artificial (2) Laminar (3) Laminar +SPS'
      read(*,*) i_viscos
      write(*,*) i_viscos

      if (i_viscos.eq.1) then
          write(*,*) 'Alpha  ??'
          read(*,*)  viscos_val
      else
          write(*,*) 'Laminar Viscosity Nu  ??'
          read(*,*)  viscos_val
      endif

      write(*,*) viscos_val


c     vorticity printing
      write(*,*) 'vorticity printing ? (1=yes)'
      read(*,*) i_vort
      write(*,*) i_vort

c     pressure terms

      print*
      write(*,*) ' Choose Equation of State (EoS)'
      write(*,*) ' 1: Weakly Compressible Fluid (Taits equation)'
      write(*,*) ' 2: Ideal Gas Equation'
      write(*,*) ' 3: Weakly Compressible Fluid (Morris equation)'
      write(*,*) ' 4: Incompressible Fluid (Poissons equation)'
      read(*,*) i_EoS
      print*,'i_EoS' ,i_EoS

c     Define pressure terms


      write(*,*) ' Maximum Depth in the simulation, h_SWL'
      read(*,*)  h_SWL
      write(*,*) h_SWL


      if(i_EoS.eq.1.or.i_EoS.eq.3) then

        write(*,*) ' Coefficient of Speed of sound ( 10 , 40 ) ??'
        read(*,*) coef
        write(*,*) coef

        coef2=coef*coef

        rho0=1000.0 !initial density
        gamma=7.0
        expont=1./gamma
        TE0=0.0

        B=coef2*rho0*g*(h_SWL)/7.0  ! B is estimated using h_SWL and coef
        Co=((B*7.0)/rho0)**.5

        write(*,*) ' B=coef2*rho0*g*(h_SWL)/7.0 =',B
        write(*,*) ' Co=',coef,'*V=',Co

      else if (i_EoS.eq.2) then

        write(*,*) ' B and coefficient of speed of sound are not used '
        coef=0.0

        rho0=1.0 
        gamma=1.4
        expont=1./gamma
        write(*,*) ' Set an intial value for TE '
        read(*,*) TE0
        write(*,*) TE0
		

      else if (i_EoS.gt.3) then

         write(*,*) 'It is not avaliable yet'
         stop

      endif

c	boundary type

      write(*,*)'Choose Boundary Type'
      write(*,*)'Monaghan & Kos (1999) = 1'
      write(*,*)'Dalrymple & Knio (2001) = 2'
      read(*,*) IBC
      write(*,*) IBC 

      if(iBC.eq.1)then
        print*,' Enter wall viscosity value for Monaghan & Kos (1999) '
        read(*,*) visc_wall
        print*,' visc_wall = ',visc_wall
        ndt_DBCPerform=0
      else
        visc_wall = 0.0
        print*,' Use Hughes & Graham (2009) correction ?  '
        print*,' ndt_DBCPerform ? (1 means no correction) '
        read(*,*) ndt_DBCPerform
        if(ndt_DBCPerform.eq.0) ndt_DBCPerform=1
        write(*,*) ndt_DBCPerform
        if(i_algorithm.eq.2.and.ndt_DBCPerform.ne.1)then
           write(*,*) 'Hughes & Graham (2009) correction is not 
     & going to work for Verlet algorithm'
           stop
        elseif(i_algorithm.eq.2) then
           write(*,*) 'Hughes & Graham (2009) correction is not 
     & avaliable for Verlet algorithm'
        endif
      endif


c       Define computational Space

      write(*,*) ' Geometry of the zone'
      write(*,*) ' (1) BOX '
      write(*,*) ' (2) BEACH '
      write(*,*) ' (3) EXTERNALLY FIXED GEOMETRY '
      read(*,*)  i_geometry
      write(*,*) i_geometry

      if (i_geometry.ne.3) then
         write(*,*) 'Fluid Particle Lattice Structure ??'
         print*,' (1) SC; (2) BCC'
         read(*,*) lattice
         write(*,*) 'Lattice Type = ', lattice
         if (lattice.eq.1) then !SC
           vnorm_mass=1.
         else  !BCC
           vnorm_mass=0.5
         endif
      endif

      call CPU_TIME (gen_Particles_time0)

      if (i_geometry.eq.1) then
         call box(g,expont,dx,dy,dz)
      else if (i_geometry.eq.2) then
         call beach(g,expont,dx,dy,dz)
      else if (i_geometry.eq.3) then
         call external_geometry(dx,dy,dz)
      else
         write(*,*) i_geometry,' is not a valid Geometry'
         stop
      endif

      call CPU_TIME (gen_Particles_time1)

      eps=0.5

      write(*,*)' please input the tmax and out:'
      read(*,*) tmax,out
      write(*,*) tmax,out
      write(*,*)' initial recording time'
      read(*,*) trec_ini
      write(*,*) trec_ini
      write(*,*)' For detailed recording during RUN'
      write(*,*)' step, start, end'
      read(*,*) dtrec_det, t_sta_det, t_end_det
      write(*,*) dtrec_det, t_sta_det, t_end_det

      time=0.
      write(*,*) 'input dt:, Variable time step (1=y)'
      read(*,*) dt,ivar_dt
      write(*,*) dt,ivar_dt

      write(*,*) 'CFL number (0.1 - 0.5)'
      read(*,*) CFL_number
      write(*,*) CFL_number

      write(*,*) ' h=coefficient*sqrt(dx*dx+dy*dy+dz*dz)'
      write(*,*) ' Coefficient ??'
      read(*,*)  coefficient
      write(*,*) coefficient




c     ______ Time Step Stability Check
      Cr = CFL_number
      Cmax = sqrt(B*gamma/rho0)  !note adjusted formula
      Vmax = max(1.1*Cmax,Cmax)
      dtmax=Cr*min(dx,dz)/Vmax

      print*
      print*,'--- Courant Stability Check ---'
      print*,'Courant No. Cr = ',Cr
      print*,'B    =           ',B
      print*,'Cmax =           ',Cmax
      print*,'dtmax: ',dtmax
      if(dtmax.lt.dt) then
        print*,'dt is TOO big!!!!'
        print*,'suggest: ',dtmax
      endif
      print*
      print*,'--- Viscous Stability Check ---'
      print*,'Kinematic Viscosity',1.0e-6
      print*,'Visc =           ',1.0e-6
      dtmax=Cr*((min(dx,dz))**2)/(101.0*1.0e-6)
      write(*,*) 'dtmax: ',dtmax
      if(dtmax.lt.dt) then
        print*,'dt is TOO big!!!!'
        print*,'suggest: ',dtmax
      endif
      print*
      write(*,*) 'dt: ',dt
      print*


c     --- Choose Riemann Solver ---
      write(*,*) ' Choose Riemann Solver'
      write(*,*) ' 0: None '
      write(*,*) ' 1: Conservative (Vila 1999)'
      write(*,*) ' 2: Non-Conservative (Parshikov et al. 2000)'
      read(*,*)  iRiemannSolver
      print*,'iRiemannSolver ' ,iRiemannSolver
      if(iRiemannSolver.eq.1)then
        i_ConservativeFormulation = 1
      else
        i_ConservativeFormulation = 0
      endif
      if(iRiemannSolver.gt.0)then
        if(i_viscos.eq.1)then
          print*,'Cannot have Riemann Solver and Artificial Viscosity'
          stop
        endif
        if(i_densityFilter.gt.0)then
          print*
          print*,'It is not advised to have a density filter'
          print*,'AND a Riemann Solver simultaneously'
          print*,'DISABLING density filter '
          i_densityFilter = 0
          print*,'i_densityFilter =',i_densityFilter
        endif
        write(*,*) ' Use TVD, slope limiter (beta_lim)?'
        read(*,*) iTVD, beta_lim
        print*,'iTVD, beta_lim ' ,iTVD, beta_lim
        Cr = CFL_number
        Cmax = sqrt(B*gamma/rho0)  !note adjusted formula
        Vmax = max(1.1*Cmax,Cmax)
        dtmax_Riemann=Cr*min(dx,dz)/Vmax
        print*
        print*,'--- Riemann Courant Stability Check ---'
        print*,'Courant No. Cr = ',Cr
        print*,'B    =           ',B
        print*,'Cmax =           ',Cmax
        print*,'dtmax_Riemann:   ',dtmax_Riemann
        if(dtmax_Riemann.lt.dt) then
          print*,'dt is TOO big!!!!'
          print*,'suggest: ',dtmax_Riemann
        endif
      else
        iTVD = 0
        beta_lim = 0.0
      end if

      dy = 0D0
      dpx = dx
      dpy = dy
      dpz = dz
188   format(es13.7)      


c	___________ File INDAT

      open(11,file='INDAT')

      write(11,*)i_kernel					! 1
      write(11,*)i_algorithm
      write(11,*)i_densityFilter
      write(11,*)i_viscos
      write(11,*)iBC						! 5
      write(11,*)i_periodicOBs(1)
      write(11,*)i_periodicOBs(2)
      write(11,*)i_periodicOBs(3)
      write(11,*)lattice
      write(11,*)i_EoS						! 10
      write(11,*)h_SWL
      write(11,*)B
      write(11,*)gamma
      write(11,*)coef
      write(11,*)eps						! 15
      write(11,*)rho0
      write(11,*)viscos_val
      write(11,*)visc_wall
      write(11,*)vlx
      write(11,*)vly						! 20
      write(11,*)vlz
      write(11,188)dx
      write(11,188)dy
      write(11,188)dz
      h=coefficient*sqrt(dx*dx+dz*dz)
      write(11,*)h						! 25
      write(11,*)np
      write(11,*)nb
      write(11,*)nbf
      write(11,*)ivar_dt
      write(11,*)dt						! 30
      write(11,*)tmax
      write(11,*)out
      write(11,*)trec_ini
      write(11,*)dtrec_det
      write(11,*)t_sta_det					! 35
      write(11,*)t_end_det
      write(11,*)i_restartRun
      write(11,*)CFL_number
      write(11,*)TE0
      write(11,*)i_kernelcorrection				! 40
      write(11,*)iRiemannSolver
	write(11,*)iTVD
	write(11,*)beta_lim
	write(11,*)i_vort
      write(11,*)ndt_VerletPerform
      write(11,*)ndt_FilterPerform
	write(11,*)ndt_DBCPerform
	close(11)
	  
      xb_min = 0.0
      xb_max = vlx
      yb_min = 0.0
      yb_max = vly
      zb_min = 0.0
      zb_max = vlz

      call position_check(dx,dz)

c	________ SUBROUTINES FOR DEALING WITH NORMALS (SERIAL Version)

      if(iBC.eq.1)then
        call normals_Calc_2D
        call normals_FileWrite_2D    
      end if
      print*

      print*,'Writing IPART file ...'
      print*

c	__________  File IPART

      open(13,file='IPART')

      do i=1,np
      
        write(13,122) xp(i),zp(i),up(i),wp(i),rhop(i),p(i),pm(i)

      enddo


122   format(7e16.8)
      close(13)

      print*,'IPART file written'
      print*

c	__________  File MATLABIN

      open(18,file='matlabin')
        write(18,*) np
        write(18,*) vlx
        write(18,*) vly		!  0 if 2D
        write(18,*) vlz
        write(18,*) out
        write(18,*) nb
        write(18,*) nbf
      close(18)


c	__________  File RESTART  
      !--  Printout RESTART file new run using CheckPointing --
      if(i_restartRun.eq.2)then
         open(44,file='RESTART',status='unknown')
         write(44,146)0,0.0,0,dt    !=itime,time,ngrab,dt
         close(44)
      endif
146   format(i8,e16.8,i8,e16.8)


c	________________________   TO COMPILE

      write(*,*) 'Which compiler is desired? '
      write(*,*) 'Linux: 1 = gfortran, 2 = ifort'
      write(*,*) 'Windows: 3 = ifort, 4 = Silverfrost FTN95'
      read(*,*) i_compile_opt
      write(*,*) i_compile_opt
      if(i_compile_opt.eq.1) then
        call tocompile_gfortran
      elseif(i_compile_opt.eq.2) then
        call tocompile_ifort
      elseif(i_compile_opt.eq.3) then
        call tocompile_win_ifort
      else  !if(i_compile_opt.eq.4) then
        call tocompile_ftn95
      end if
      
c     ______  Choose Precision of XYZ Variables for Run _________
      call precisionWrite(i_compile_opt)

      print*,'nbf  ',nbf
      print*,'nbfm ',nbfm
      print*,'nb   ',nb
      print*,'np   ',np

      write(*,145)'Generated np = ',np,' particles in ',
     &        gen_Particles_time1-gen_Particles_time0,' seconds'
145   format(a15,1x,i9,a14,1x,f5.2,a8)

      stop

      end    ! END OF PROGRAM


c	________________________   SUBROUTINE BOX
c
      subroutine box(g,expont,dx,dy,dz)
      include 'common.gen2D'


       pi=4.*atan(1.)
c       
c     _______ Generate a 2D BOX
c
    
 
         print*
         print*,'--- Generating 2-D Box ---'
         print*
 
         write(*,*) ' Box dimension LX,LZ?'		
         read(*,*) vlx,vlz
         write(*,*) vlx,vlz
         vly=0.
       
         write(*,*) ' Spacing dx,dz?'			
         read(*,*)  dx,dz
         write(*,*) dx,dz
         dy=0.
       
         write(*,*) ' Inclination of floor in X ( beta ) ?? '  
         read(*,*)  beta_deg
         write(*,*) beta_deg
         beta=beta_deg*pi/180
       
         theta=0.
       
         N=nint(vlx/dx)
         M=0		
         L=nint(vlz/dz)

       
       write(*,*) ' Dimension (pixels)', N,M,L
       
       print*
       write(*,*) ' Enter Periodicity Information:'
       write(*,*) ' Add Periodic Lateral boundaries'
       write(*,*) ' in X, Y, & Z-Directions ? (1=yes)'
       read(*,*) i_periodicOBs(1),i_periodicOBs(2),i_periodicOBs(3)
       print*
       print*,'i_periodicOBs'
       print*,'X-Direction: ',i_periodicOBs(1)
       print*,'Y-Direction: ',i_periodicOBs(2)
       print*,'Z-Direction: ',i_periodicOBs(3)
       if((i_periodicOBs(1).ne.0.and.i_periodicOBs(1).ne.1).or.
     &    (i_periodicOBs(2).ne.0.and.i_periodicOBs(2).ne.1).or.
     &    (i_periodicOBs(3).ne.0.and.i_periodicOBs(3).ne.1))then    
          print*
          print*,'Not a valid Periodic Boundaries Option'
          print*
          stop
       end if

      nn=0	! Initial number of particles

      
       if(i_periodicOBs(1).ne.1)then  !Lateral Boundaries
         call boundaries_left(nn,0,N1,M,L,dx,dy,dz,beta,theta)     !x = x1
         call boundaries_right(nn,0,N,M,L,dx,dy,dz,beta,theta)     !x = x2
       endif
       
       if(i_periodicOBs(3).ne.1)then  !Lateral Boundaries
         call boundaries_bottom(nn,0,N,M,L,dx,dy,dz,beta,theta)    !z = z1
       endif


	write(*,*) ' Add wall inside the box (1=y)'
	read(*,*)  iopt_addwall
	write(*,*) iopt_addwall

	do while (iopt_addwall.eq.1)
	   call wall(nn,dx,dy,dz,beta,theta)
	   write(*,*) ' Add another wall inside the box (1=y)'
	   read(*,*)  iopt_addwall
	   write(*,*) iopt_addwall
	enddo


      open(55,file='obstacle')

	write(*,*) ' Add obstacle inside the box (1=y)'
	read(*,*)  iopt_obst
	write(*,*) iopt_obst

      write(55,*) iopt_obst

      n_trapezoid = 0
      do while (iopt_obst.eq.1)
          i_obstacleShape = 0
          write(*,*) 'Choose rectangular (1) or trapezoid (2) '
          read(*,*) i_obstacleShape
          write(*,*)'i_obstacleShape ',i_obstacleShape
	   if(i_obstacleShape.eq.1)then
            call obstacle (nn,dx,dy,dz,beta,0.,iopt_obst)
          elseif(i_obstacleShape.eq.2)then
            call trapezoid(nn,dx,dy,dz,beta,0.,iopt_obst)
          else
            print*,'Invalid option'
            stop
          endif
          write(*,*) ' Add another obstacle inside the beach (1=y)'
          read(*,*)  iopt_obst
          write(*,*) iopt_obst
      enddo

	close(55)

	nbf=nn

	write(*,*) 'No of B. Fixed Particles',nbf


	open(66,file='wavemaker')

	write(*,*) ' Add wavemaker inside the box (1=y)'
	read(*,*)  iopt_wavemaker
	write(*,*) iopt_wavemaker

	write(66,*) iopt_wavemaker

	do while (iopt_wavemaker.eq.1)
	   write(*,*) ' Just in Y direction'
	   call wavemaker(nn,dx,dy,dz,beta,theta)
	   write(*,*) ' Add another wavemaker inside the box (1=y)'
	   read(*,*)  iopt_wavemaker
	   write(*,*) iopt_wavemaker
	enddo

	close(66)


      open(77,file='gate')

	write(*,*) ' Add gate inside the box (1=y)'
	read(*,*)  iopt_gate
	write(*,*) iopt_gate

	write(77,*) iopt_gate

	do while (iopt_gate.eq.1)
	   write(*,*) ' Just in Y direction'
	   call gate(nn,dx,dy,dz,beta,theta)
	   write(*,*) ' Add another gate inside the box (1=y)'
	   read(*,*)  iopt_gate
	   write(*,*) iopt_gate
	enddo

	close(77)
       open(88,file='Tsunami_Landslide.txt')
       iopt_RaichlenWedge=0
       write(88,*)iopt_RaichlenWedge
       if(iopt_RaichlenWedge.eq.1.and.iopt_wavemaker.eq.1)then
         print*,'Cannot have Raichlen Wedge and wavemaker'
         print*,' at the same time'
         stop
       end if
       if(iopt_RaichlenWedge.eq.1)then
          call RaichlenWedge_Particles(nn,N1,N,M,dx,dy,dz,beta)
       endif
       close(88)

       print*
       num_FB = 0
       open(88,file='Floating_Bodies.txt')
       write(*,*)' Add Floating Bodies (1=yes) ?'
       read(*,*)  iopt_FloatingBodies
       write(*,*) iopt_FloatingBodies
       write(88,*)iopt_FloatingBodies
       nbfm = nn
       print*,'nbfm ',nbfm
       write(88,*)nbfm
       do while (iopt_FloatingBodies.eq.1)
         call FloatingBodies_Particles(nn,N1,N,M,dx,dy,dz,beta)
         write(*,*) ' Add another Floating Body (1=yes)'
         read(*,*)  iopt_FloatingBodies
         write(*,*) iopt_FloatingBodies
       enddo
       !write(88,*)num_FB
       close(88)

	nb=nn

	write(*,*) 'No of B. Particles',nb

	call fluid_particles(nn,N,M,L,dx,dy,dz,expont,g,i_correct_pb,
     +            XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax,beta,theta)

	np=nn
	write(*,*) 'No of Particles',np

	if(i_correct_pb.eq.1) then
	   call correct_p_boundaries(expont,g,XXmin,XXmax,YYmin,
     +			             YYmax,ZZmin,ZZmax,dx,dy,dz)
        else
	   call p_boundaries(dx,dy,dz)
	endif

	write(*,*) 'End'
	end


c	________________________   SUBROUTINE BEACH

      subroutine beach(g,expont,dx,dy,dz)
      include 'common.gen2D'


	pi=4.*atan(1.)

             
c
c      |                                                         /
c      |                                                        /
c      |                                                       /
c      |                                                      /
c      |                                                     /
c      |                                                    /
c      |                                                   /
c      |_________________________________________________ /
c      |                                                 /
c      |                                                /
c      |                                               /
c      |          i_water = 1                         /
c      |                                             /
c      |                                            /
c      |                                           /
c      |                                          /    angle beta
c      |_________________________________________/ _ _ _ _ _
c
c      <--------------- vlx1 = N1.dx------------><-- vlx2 = N2.dx -->
c      <------------------------- vlx = N.dx------------------------>
c


c     _________ Generate a BEACH

   

         print*
         print*,'--- Generating 2-D Beach ---'
         print*
         write(*,*) ' Box dimension LX,LZ?'   
         read(*,*)  vlx,vlz
         write(*,*) vlx,vlz
       
         write(*,*) ' Spacing dx,dz?'
         read(*,*)  dx,dz
         write(*,*) dx,dz
         dy = 0.
       
         write(*,*) ' Length of the flat bed before beach ?? '  
         read(*,*)  vlx1
         write(*,*) vlx1
       
         vlx2=vlx-vlx1
       
         print*,'vlx1 ',vlx1
         print*,'vlx2 ',vlx2
         print*,'vlx  ',vlx
       
         write(*,*) ' Slope (deg) of the inclined plane ( beta ) ?? '  
         read(*,*)  beta_deg
         write(*,*) 'beta_deg ',beta_deg
         beta = beta_deg*pi/180.0
       
         N =nint(vlx/dx)
         N1=nint(vlx1/dx)
         N2=nint(vlx2/dx)
         M=0
         L =nint(vlz/dz)

         print*
         write(*,*) ' Enter Periodicity Information:'
         write(*,*) ' Add Periodic Lateral boundaries'
         write(*,*) ' in X, Y, & Z-Directions ? (1=yes)'
         read(*,*) i_periodicOBs(1),i_periodicOBs(2),i_periodicOBs(3)
         print*
         print*,'i_periodicOBs'
         print*,'X-Direction: ',i_periodicOBs(1)
         print*,'Y-Direction: ',i_periodicOBs(2)
         print*,'Z-Direction: ',i_periodicOBs(3)
         if(i_periodicOBs(1).eq.1.or.i_periodicOBs(3).eq.1)then
           print*
           print*,'Periodic Boundaries in X & Z-Directions'
           print*,'Not yet activated in code'
           print*
           stop
         else if((i_periodicOBs(1).ne.0.and.i_periodicOBs(1).ne.1).or.
     &           (i_periodicOBs(2).ne.0.and.i_periodicOBs(2).ne.1).or.
     &           (i_periodicOBs(3).ne.0.and.i_periodicOBs(3).ne.1))then    
           print*
           print*,'Not a valid Periodic Boundaries Option'
           print*
           stop
         end if
       

       
      write(*,*) ' Dimension (pixels)', N,M,L
       
	nn=0	! Initial number of particles


       write(*,*)' If wavemaker will be added, left wall is not needed'
       read(*,*)  iopt_wavemaker
       write(*,*) iopt_wavemaker


c     ______ plane
                       
       if(i_periodicOBs(1).ne.1)then  
         if(iopt_wavemaker.ne.1)then
           call boundaries_left(nn,0,N1,M,L,dx,dy,dz,0.,0.)      !x = x1
         endif
       endif

c     _______ tilted plane
       

       if(i_periodicOBs(1).ne.1)then  
         call boundaries_right (nn,N1,N,M,L,dx,dy,dz,beta,0.)    !x = x2
       endif
       if(i_periodicOBs(3).ne.1)then  !Lateral Boundaries
         call boundaries_bottom(nn,0,N1,M,L,dx,dy,dz,0.  ,0.)    !flat bottom plane
         if(beta.gt.0.)then
           call boundaries_bottom(nn,N1,N,M,L,dx,dy,dz,beta,0.)    !inclined beach plane
         endif
       endif
     

      open(55,file='obstacle')

      write(*,*) ' Add obstacle inside the beach (1=y)'
      read(*,*)  iopt_obst
      write(*,*) iopt_obst

      write(55,*) iopt_obst

      n_trapezoid = 0
      do while (iopt_obst.eq.1)
          i_obstacleShape = 0
          write(*,*) 'Choose rectangular (1) or trapezoid (2) '
          read(*,*) i_obstacleShape
          write(*,*)'i_obstacleShape ',i_obstacleShape
	   if(i_obstacleShape.eq.1)then
            call obstacle (nn,dx,dy,dz,beta,0.,iopt_obst)
          elseif(i_obstacleShape.eq.2)then
            call trapezoid(nn,dx,dy,dz,beta,0.,iopt_obst)
          else
            print*,'Invalid option'
            stop
          endif
          write(*,*) ' Add another obstacle inside the beach (1=y)'
          read(*,*)  iopt_obst
          write(*,*) iopt_obst
      enddo

	close(55)

	nbf=nn

	write(*,*) 'No of B. Fixed Particles',nbf

	open(66,file='wavemaker')

	write(66,*) iopt_wavemaker

      iopt_wavemaker_temp = iopt_wavemaker  !Using temporary variable as iopt_wavemaker is used in fill_part
      do while(iopt_wavemaker_temp.eq.1)
         write(*,*) ' Just in Y direction'
         call wavemaker(nn,dx,dy,dz,beta,0.)
         write(*,*) ' Add another wavemaker For beach (1=y)'
         read(*,*)  iopt_wavemaker_temp
         write(*,*) iopt_wavemaker_temp
      enddo

	close(66)

      open(77,file='gate')

	write(*,*) ' Add gate inside the beach (1=y)'
	read(*,*)  iopt_gate
	write(*,*) iopt_gate

	write(77,*) iopt_gate

	do while (iopt_gate.eq.1)
	   write(*,*) ' Just in Y direction'
	   call gate(nn,dx,dy,dz,beta,0.)
	   write(*,*) ' Add another gate inside the beach (1=y)'
	   read(*,*)  iopt_gate
	   write(*,*) iopt_gate
	enddo

	close(77)
       
      open(88,file='Tsunami_Landslide.txt')
       write(*,*)' Add slidingWedge ?'
       read(*,*)  iopt_RaichlenWedge
       write(*,*) iopt_RaichlenWedge
       write(88,*)iopt_RaichlenWedge
       if(iopt_RaichlenWedge.eq.1.and.iopt_wavemaker.eq.1)then
         print*,'Cannot have Raichlen Wedge and wavemaker'
         print*,' at the same time'
         stop
       end if
       if(iopt_RaichlenWedge.eq.1)then
          call RaichlenWedge_Particles(nn,N1,N,M,dx,dy,dz,beta)
       endif
       close(88)

       print*
       num_FB = 0
       open(88,file='Floating_Bodies.txt')
       write(*,*)' Add Floating Bodies (1=yes) ?'
       read(*,*)  iopt_FloatingBodies
       write(*,*) iopt_FloatingBodies
       write(88,*)iopt_FloatingBodies
       nbfm = nn
       print*,'nbfm ',nbfm
       write(88,*)nbfm
       do while (iopt_FloatingBodies.eq.1)
         call FloatingBodies_Particles(nn,N1,N,M,dx,dy,dz,beta)
         write(*,*) ' Add another Floating Body (1=yes)'
         read(*,*)  iopt_FloatingBodies
         write(*,*) iopt_FloatingBodies
       enddo
       !write(88,*)num_FB
       close(88)
       
       nb=nn
       write(*,*) 'No of B. Particles',nb

      i_correct_pb=1

c	______ water in the flat region

      write(*,*) ' Add water in the flat region ?? (1=y)'
      read(*,*)    i_water
      write(*,*)   i_water

      if(i_water.eq.1) then
         call fill_part(nn,XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax,
     +               g,dx,dy,dz,expont,0.,0.)
      endif

      i_correct_pb=1

c     _______ water in the inclined region

      write(*,*) ' Add water in the inclined region ?? (1=y)'
      read(*,*)    i_water2
      write(*,*)   i_water2

      if(i_water2.eq.1) then
         call fill_part(nn,XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax,
     +               g,dx,dy,dz,expont,beta,0.)
      endif


c	profile of the wave

      write(*,*) ' Add a solitary wave ?? (1=y)'
      read(*,*)    i_wave
      write(*,*)   i_wave

      if (i_wave.eq.1) then
          call wave(nn,XXmin,XXmax,YYmin,YYmax
     +                ,ZZmin,ZZmax,g,dx,dy,dz,expont)
      endif

      np=nn
      write(*,*) 'No of Particles',np


      if(i_correct_pb.eq.1) then
         call correct_p_boundaries(expont,g,XXmin,XXmax,YYmin,
     +                             YYmax,ZZmin,ZZmax,dx,dy,dz)
      else
         call p_boundaries(dx,dy,dz)
      endif

      return
      end

c	_________________ SUBROUTINE EXTERNAL_GEOMETRY

	subroutine external_geometry(dx,dy,dz)
      include 'common.gen2D'


      i_correct_pb=1
	

	    open(10,file='INI_POS_VEL')
	    read(10,*) np,nb,nbf
	    read(10,*) dx,dz
	    read(10,*) lattice ! (1) SC; (2) BCC

	    if (lattice.eq.1) then !SC
	       vnorm_mass=1.
	    else  !BCC
	       vnorm_mass=0.5
	    endif

	    do n=1,np
	       read(10,*) xp(n),zp(n),up(n),wp(n)
	       call pressure(n,zp(n),dx,0.,dz,expont,g)
	    enddo
	    close(10)

	    XXmax=-1.e+10
	    ZZmax=-1.e+10

	    XXmin=1.e+10
	    ZZmin=1.e+10

	   do im=1,nb ! Recalculate the extent of the medium
	       XXmax=max(XXmax,xp(im))
	       XXmin=min(XXmin,xp(im))
	       ZZmax=max(ZZmax,zp(im))
	       ZZmin=min(ZZmin,zp(im))
	   enddo
	   if(i_correct_pb.eq.1) then
	      call correct_p_boundaries(expont,g,XXmin,XXmax,0.,
     +                                0.,ZZmin,ZZmax,dx,0.,dz)
	   else
	       call p_boundaries(dx,dy,dz)
	   endif

 
      return
	end



c     _________________ SUBROUTINE BOUNDARIES LEFT

      subroutine boundaries_left(NN,N0,N,M,L,dx,dy,dz,beta,theta)
      include 'common.gen2D'

      ! Lateral Walls: Left  x = x1
      print*,'Left  wall (x = x1) starts at ',nn+1
      N_start  = N0
      N_finish = N
      L_finish = L

      if(iBC.eq.2)then        !Dalrymple BC - Fixed Fluid Particles

        !- Inner Layer -
        i=N0
        do j=0,M
          hk=(i-N0)*dx*tan(beta)+j*dy*tan(theta)
          Lini=nint(hk/dz)
          do k=Lini,L
            nn=nn+1
            call pos_veloc(nn,i*dx,j*dy,k*dz,0.,0.,0.)
          enddo
        enddo

	 !- Outer Layer -
 
          do j=0,M
            i=N0
            hk=(i-N0)*dx*tan(beta)+j*dy*tan(theta)
            Lini=nint(hk/dz)
            do k=Lini,L
              nn=nn+1
              call pos_veloc(nn,(i-0.5)*dx,j*dy,(k-0.5)*dz,0.,0.,0.)
            enddo
          enddo
  
        
      else if(iBC.eq.1)then   !Repulsive BC - Repulsive Force Particles with Normals
        
        if(i_periodicOBs(2).eq.1)then
          M_start = 0 
          M_finish = M
          
        else
          M_start = 0 
          M_finish = M
        end if
        N_start = 0
        N_finish = 0
        !- One Layer -
        i=N_start
        do j=M_start,M_finish    !In 2-D both are zero --> one column of BPs
          hk=(i-N_start)*dx*tan(beta)+j*dy*tan(theta)
          L_start=nint(hk/dz)
          do k=L_start,L_finish
            nn=nn+1
            call pos_veloc(nn,i*dx,j*dy,k*dz,0.,0.,0.)
           !2-D Normals
              if(k.eq.L_start)then
                i_minus1 = nn + 1
                i_plus1  = nn 
              else if(k.eq.L_finish)then
                i_minus1 = nn 
                i_plus1  = nn - 1
              else
                i_minus1 = nn + 1
                i_plus1  = nn - 1
              end if
              !-- Normal information - neighbour data for        --
              !-- Repulsive Boundary Particles (BPs)             --
              iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
              iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          
          enddo
        enddo

      end if     !End of:  if(iBC.eq.1)then

      return
      end

c     _________________ SUBROUTINE BOUNDARIES RIGHT

      subroutine boundaries_right(NN,N0,N,M,L,dx,dy,dz,beta,theta)
      include 'common.gen2D'

      ! Lateral Walls: Right  x = x2
      print*,'Right wall (x = x2) starts at ',nn+1
      
      N_start  = N0
      N_finish = N
      L_finish = L
      if(i_periodicOBs(2).eq.1)then
        M_start = 0
        M_finish = M
      else
        M_start = 0 
        M_finish = M
      end if

      if(iBC.eq.2)then        !Dalrymple BC - Fixed Fluid Particles

        !- Inner Layer -
        do j=M_start,M_finish
          i=N
          hk=(i-N0)*dx*tan(beta)+j*dy*tan(theta)
          Lini=nint(hk/dz)
          do k=Lini,L
            nn=nn+1
            call pos_veloc(nn,i*dx,j*dy,k*dz,0.,0.,0.)
          enddo
        enddo

	 !- Outer Layer -
        
          do j=M_start,M_finish
            i=N
            hk=(i-N0)*dx*tan(beta)+j*dy*tan(theta)
            Lini=nint(hk/dz)
            do k=Lini,L
              nn=nn+1
              call pos_veloc(nn,(i+0.5)*dx,j*dy,(k-0.5)*dz,0.,0.,0.)
            enddo
          enddo
  
        
      else if(iBC.eq.1)then   !Repulsive BC - Repulsive Force Particles with Normals
        
        N_start = N0
        N_finish = N
        !- One Layer -
        i=N_finish
        do j=M_start,M_finish      !In 2-D both are zero --> one column of BPs
          hk=(i-N_start)*dx*tan(beta)+j*dy*tan(theta)
          L_start=nint(hk/dz)
          do k=L_start,L_finish
            nn=nn+1
            call pos_veloc(nn,i*dx,j*dy,k*dz,0.,0.,0.)
    
              if(k.eq.L_start)then
                i_minus1 = nn
                i_plus1  = nn + 1 
              else if(k.eq.L_finish)then
                i_minus1 = nn - 1 
                i_plus1  = nn
              else
                i_minus1 = nn - 1
                i_plus1  = nn + 1
              end if
              !-- Normal information - neighbour data for        --
              !-- Repulsive Boundary Particles (BPs)             --
              iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
              iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
    
          enddo
        enddo
        
      end if     !End of:  if(iBC.eq.1)then

      return
      end


c     _________________ SUBROUTINE BOUNDARIES BOTTOM

      subroutine boundaries_bottom(NN,N0,N,M,L,dx,dy,dz,beta,theta)
      include 'common.gen2D'

      double precision hk

      !- Bottom -
      if(beta.gt.0.0)then
        print*,'Inclined Plane boundary starts at ',nn+1
      else 
        print*,'Bottom boundary starts at         ',nn+1
      end if

      !- One Layer -
      if(i_periodicOBs(2).eq.1)then
        M_start = 0 
        M_finish = M
       
      else
        M_start = 0 
        M_finish = M
      end if
      
      
      if(iBC.eq.2)then        !Dalrymple BC - Fixed Fluid Particles

        N_start = N0
        N_finish = N-1

        if(i_periodicOBs(1).eq.1)then
          !N_start = N_start-1
          N_finish = N_finish+1
        endif

        !- Inner Layer -
        do i=N_start+1,N_finish   !N0+1,N-1     !do i=N0,N   <--- Be careful here not to place particles ontop of each other
          do j=M_start,M_finish
            hk=(i-N0)*dx*tan(beta)+j*dy*tan(theta)
            nn=nn+1
            call pos_veloc(nn,i*dx,j*dy,hk,0.,0.,0.)
          enddo
        enddo

        if(i_periodicOBs(1).eq.1)then
          !N_start = N_start+1
          N_finish = N_finish-1
        endif

        !- Outer Layer -

          do i=N_start,N_finish   !N0,N-1
            do j=0,M
              hk=(i-N0)*dx*tan(beta)+j*dy*tan(theta)
              nn=nn+1
              call pos_veloc(nn,(i+0.5)*dx,j*dy,hk-0.5*dz,0.,0.,0.)
            enddo
          enddo

        
      else if(iBC.eq.1)then   !Repulsive BC - Repulsive Force Particles with Normals
        

        if(i_geometry.eq.2)then   !Sloping Beach
          !N_finish = int((x_beach_start - 0.0)/dx) + 1
          if(beta.gt.0.0)then  !Inclined Plane Section
            N_start = int((vlx1/dx)) + 1
            N_finish = N
          else                     !Flat Bed Section
            N_start = 0
            N_finish = int((vlx1/dx))
          end if
        else
          N_start = 0
          N_finish = N
        end if

        do i=N_start,N_finish
          do j=M_start,M_finish         !do j=0,M

            hk=(i-N_start)*dx*tan(beta)+j*dy*tan(theta)  
            nn=nn+1
            call pos_veloc(nn,i*dx,j*dy,hk,0.,0.,0.)

              if(i.eq.N_start)then
                i_minus1 = nn
                i_plus1  = nn + 1
              else if(i.eq.N_finish)then
                i_minus1 = nn - 1
                i_plus1  = nn
              else
                i_minus1 = nn - 1
                i_plus1  = nn + 1
              end if
              !-- Normal information - neighbour data for        --
              !-- Repulsive Boundary Particles (BPs)             --
              iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
              iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
  
          enddo
        enddo
        
      end if     !End of:  if(iBC.eq.1)then

      return
      end

c      ____

c	___________________   SUBROUTINE POS_VELOC

       subroutine pos_veloc(nn,xpos,ypos,zpos,uveloc,vveloc,wveloc)
       include 'common.gen2D'

       double precision xpos,zpos
	
	xp(nn)=xpos
c	yp(nn)=ypos
	zp(nn)=zpos
	up(nn)=uveloc
c	vp(nn)=vveloc
	wp(nn)=wveloc

	return
	end

C	____

c	____________________ SUBROUTINE FLUID_PARTICLES

	subroutine fluid_particles(nn,N,M,L,dx,dy,dz,expont,g,
     +                              i_correct_pb,XXmin,XXmax,YYmin,
     +                              YYmax,ZZmin,ZZmax,beta,theta)
      include 'common.gen2D'
      
	write(*,*) 'Initial conditions'
	write(*,*) '1) Set of particles without grid'
	write(*,*) '2) Fill a square shaped region '
	write(*,*) '3) Round shaped drop  '
	write(*,*) '  Which Initial condition ??'
	read(*,*)  initial
	write(*,*) initial
	 
C	 I_CORRECT_PB =   1 bound. pressure will  be corrected
C		          0 boundary pressure will not be corrected


	if(initial.eq.1) then
          i_correct_pb=0
	    call set(nn,g,expont,ZZmin,ZZmax,dx,dy,dz)
	else if(initial.eq.2) then
          write(*,*) ' Correct pressure at boundaries ?? (1=y)'
	    read(*,*)    i_correct_pb
	    write(*,*)   i_correct_pb
	    iopt_fill_part=1
	    do while (iopt_fill_part.eq.1)
	       call fill_part(nn,XXmin,XXmax,YYmin,YYmax
     +                     ,ZZmin,ZZmax,g,dx,dy,dz,expont,beta,theta)
	       write(*,*)'Last Particle',nn
             write(*,*) ' Fill a new region ?? (1=y)'
	       read(*,*)    iopt_fill_part
	       write(*,*)   iopt_fill_part
	    enddo
	else if(initial.eq.3) then
          i_correct_pb=0
	    call drop(nn,g,dx,dy,dz,expont)
	
	else
          write(*,*) initial,' is not a valid I. Condition'
	    stop
	endif

	return

	end

c	_____

c	____________________ SUBROUTINE PRESSURE

      subroutine pressure(nn,ZZmax,dx,dy,dz,expont,g)
      include 'common.gen2D'

     
      rhop(nn)=rho0*(1.0+rho0*g*(ZZmax-zp(nn))/B)**expont
      p(nn)=B*((rhop(nn)/rho0)**gamma-1.0)
	
	pm(nn)=vnorm_mass*rhop(nn)*dx*dz

      return
	end


c	____________________ SUBROUTINE CORRECT_P_BOUNDARIES

      subroutine correct_p_boundaries (expont,g,XXmin,XXmax,YYmin,
     +			                 YYmax,ZZmin,ZZmax,dx,dy,dz)
      include 'common.gen2D'

	dr=dx*dz


	do nn=1,nb
	   if (xp(nn).le.XXmax. and. xp(nn).ge.XXmin.and.
     +	       Zp(nn).le.ZZmax. and. zp(nn).ge.ZZmin) then

               rhop(nn)=rho0*(1.0+rho0*g*(ZZmax-zp(nn))/B)**expont
	   else
		     rhop(nn)=rho0
	   endif

         p(nn)=B*((rhop(nn)/rho0)**gamma-1.0)

	   if (IBC.eq.1) then
            pm(nn)=rhop(nn)*dr
	   else if (IBC.eq.2) then
            if(nn.gt.nbfm) then
               pm(nn)=rhop(nn)*dr ! only one layer of particles for Floating Bodies faces
            else
               pm(nn)=0.5*rhop(nn)*dr
            endif
	   endif

	enddo

      return

	end

c	______

c	____________________ SUBROUTINE P_BOUNDARIES

      subroutine p_boundaries (dx,dy,dz)
      include 'common.gen2D'

	dr=dx*dz


	do nn=1,nb
	   rhop(nn)=rho0
         p(nn)=0.

	   if (IBC.eq.1) then
            pm(nn)=rhop(nn)*dr
	   else if (IBC.eq.2) then
            if(nn.gt.nbfm) then
               pm(nn)=rhop(nn)*dr   ! only one layer of particles for Floating Bodies faces
            else
               pm(nn)=0.5*rhop(nn)*dr
            endif
	   endif

	enddo

      return

	end

c	______

c	____________________ SUBROUTINE SET

      subroutine set (nn,g,expont,ZZmin,ZZmax,dx,dy,dz)
      include 'common.gen2D'

      
      write(*,*) ' Number of particles ??'
	read(*,*) number
	write(*,*) number

	do num=1,number
	   nn=nn+1

	   write(*,*) 'Position X,Z'     
	   read(*,*)  xp(nn),zp(nn)
	   write(*,*) xp(nn),zp(nn)  
	   write(*,*) 'Velocity X,Z'     
	   read(*,*)  up(nn),wp(nn)
	   write(*,*) up(nn),wp(nn)	   

         call pressure(nn,zp(nn),dx,dy,dz,expont,g)
	enddo

c	NOTE:For a set of particles without grid B is estimated
c          using the Z coordinate of the box

      return
	end
c	______
c	____________________ SUBROUTINE FILL_PART
       
       subroutine fill_part(nn,XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax,
     +	                     g,dx,dy,dz,expont,beta,theta)
       include 'common.gen2D'
       

      double precision x1,z1,XXmin_initial


       
       write(*,*) ' Cube containing particles'
       write(*,*) ' XMin, Xmax ??'
       read(*,*) XXmin,XXmax
       write(*,*) XXmin,XXmax
       
       YYmax=0.
       YYmin=0.
        
       write(*,*) ' ZMin, Zmax ??'
       read(*,*) ZZmin,ZZmax
       write(*,*) ZZmin,ZZmax
       
       N_end=nint(XXmax/dx)
       N_ini=nint(XXmin/dx)
       
       Nstepsx=nint((XXmax-XXmin)/dx)+1
       Nstepsy=1
       Nstepsz=nint((ZZmax-ZZmin)/dz)+1
       
       !- Info to stop WPs being placed at same exact location as BPs -
       XXmin_initial = XXmin
       XXmax_initial = XXmax
       ZZmin_initial = ZZmin
       ZZmax_initial = ZZmax
       write(*,*)
       write(*,*)'XXmin_initial ',XXmin_initial
       write(*,*)'XXmax_initial ',XXmax_initial
       write(*,*)'ZZmin_initial ',ZZmin_initial
       write(*,*)'ZZmax_initial ',ZZmax_initial
       write(*,*)
       !---------------------
       
       if(num_FB.gt.0)then
         print*,'num_FB ',num_FB
         print*,'x_cyl_min(num_FB) ',x_cyl_min(num_FB)
         print*,'z_cyl_min(num_FB) ',z_cyl_min(num_FB)
         print*,'x_cyl_max(num_FB) ',x_cyl_max(num_FB)
         print*,'z_cyl_max(num_FB) ',z_cyl_max(num_FB)
       endif
       
       XXmin = N_ini*dx
       XXmax = N_end*dx
       
       if(XXmin.gt.XXmin_initial + vnorm_mass*0.5*dx)then
         XXmin = XXmin_initial   !Undoes rounding effect nint above for certain domains
       else if(XXmin.lt.XXmin_initial + 0.99*vnorm_mass*0.5*dx)then
         XXmin = XXmin_initial   !Undoes rounding effect nint above for certain domains
       end if
                  
       M_end=0
       M_ini=0
       
       YYmin=0.
       YYmax=0.
            
       L_end=nint(ZZmax/dz)
       ZZmax=L_end*dz
       L_ini=nint(ZZmin/dz)
       ZZmin=L_ini*dz
       
       if(ZZmax.gt.h_SWL) ZZmax = h_SWL
       
       if(iBC.eq.1)then    !Repulsive force Boundary Condition             
         L_ini = L_ini-1
         ZZmin = (L_ini+0.5)*dz
       end if
                  
       write(*,*) 'Recalculated max and min X ',XXmax,XXmin
       write(*,*) 'Recalculated max and min Y ',YYmax,YYmin
       write(*,*) 'Recalculated max and min Z ',ZZmax,ZZmin
       
       print*
       print*,'i_obstacleShape ',i_obstacleShape
       if(i_obstacleShape.gt.1)then
          print*,'n_trapezoid ',n_trapezoid
c         print*,'x,z_trapezoid(n_trapezoid,1) ',
c     &   x_trapezoid(n_trapezoid,1), z_trapezoid(n_trapezoid,1)
c         print*,'x,z_trapezoid(n_trapezoid,2) ',
c     &   x_trapezoid(n_trapezoid,2), z_trapezoid(n_trapezoid,2)
c         print*,'x,z_trapezoid(n_trapezoid,3) ',
c     &   x_trapezoid(n_trapezoid,3), z_trapezoid(n_trapezoid,3)
c         print*,'x,z_trapezoid(n_trapezoid,4) ',
c     &   x_trapezoid(n_trapezoid,4), z_trapezoid(n_trapezoid,4)
       endif
       print*,'iopt_wavemaker  ',iopt_wavemaker
       print*,'i_paddleType    ',i_paddleType
       print*,'beta',beta  
       print*    
	    
       write(*,*) ' Npnts x/y/z ', nstepsx,nstepsy,nstepsz
       if(beta.gt.0)L_end = L_ini + 10  
       
       do i=N_ini,N_end
         do j=1,Nstepsy
         
           hk=(i-N_ini)*dx*tan(beta)+j*dy*tan(theta)   
           
           if(i.eq.N_ini.and.hk.lt.ZZmin)then
             ZZmin_BC1 = hk
           end if
           
           !print*,'i,x1 ',i,XXmin+(i-N_ini)*dx           
    
c           if(iBC.eq.2)then
c             if(beta.gt.0.0)hk = hk + 0.5*dz    !Stops water particles being placed ontop of sloping bed particles
c             ZZinf=max(ZZmin_BC1,hk)
c             L_ini=nint(ZZinf/dz)
c             L_end=nint(ZZmax/dz)
c           else if(iBC.eq.1)then    !Repulsive force Boundary Condition
             ZZinf=min(ZZmin_BC1,hk)  !MIN ensures continuity of initial particle structure
             L_ini=int(ZZinf/dz)
             L_end=int(ZZmax/dz)
             if(L_end*dz.lt.ZZmax_initial)then
               L_end = L_end + 1         
             end if
c           end if
           
            
           do k=L_ini,L_end
              
             x1=XXmin+(i-N_ini)*dx
             y1=YYmin+(j-1)*dy
             y1=0
             
             z1 = ZZinf+(k-L_ini)*dz
       
             !Check limits for objects
             if(iopt_wavemaker.eq.1)then
               if(i_paddleType.eq.2)then
                 XX_extra =  X_PaddleCentre 
     &               +0.5*(1.0+(z1 - paddle_SWL)/
     &                    (flap_length+paddle_SWL))*stroke
               elseif(i_paddleType.eq.3)then
                 XX_extra =  X_PaddleStart_3
               else
                 XX_extra =  X_PaddleCentre
               end if
               if(x1.gt.XX_extra + 0.5*dx)then
                 i_createParticle = 1
               else
                 i_createParticle = 0
               end if
             else if(iopt_RaichlenWedge.eq.1)then
               if(    x1.ge.block_xstart_point-0.25*dx
     &           .and.x1.le.block_xfinish_point+0.25*dx
     &           .and.z1.le.block_top+0.5*dz)then
                 i_createParticle = 0
               else
                 i_createParticle = 1
               end if
             else
               i_createParticle = 1
             end if
             z_clearanceFactor = 0.75  !YOU CAN CHANGE THIS!!
             if(i_obstacleShape.gt.1.and.i_createParticle.eq.1)then
              i_trapezoid = 0
c              do while(i_createParticle.eq.1
c     &             .or.i_trapezoid.lt.n_trapezoid)
              do while(i_trapezoid.lt.n_trapezoid)
               i_trapezoid = i_trapezoid + 1
               if(x1.ge.x_trapezoid(i_trapezoid,4))then
                 i_createParticle = 1
               elseif(x1.ge.x_trapezoid(i_trapezoid,3))then
                 trapezoid_slope1 =
     &           (z_trapezoid(i_trapezoid,4)-z_trapezoid(i_trapezoid,3))
     &          /(x_trapezoid(i_trapezoid,4)-x_trapezoid(i_trapezoid,3))            
                 z_trap_temp = z_trapezoid(i_trapezoid,3) + 
     &           (x1 - x_trapezoid(i_trapezoid,3))*trapezoid_slope1
                 !print*,'
                 if(z1.lt.z_trap_temp+(z_clearanceFactor+0.5)*dz)then
                   i_createParticle = 0
                   i_trapezoid = n_trapezoid + 1
                 else
                   i_createParticle = 1
                 endif
               elseif(x1.ge.x_trapezoid(i_trapezoid,2))then
                 z_trap_temp =
     &           (int(z_trapezoid(i_trapezoid,2)/dz) + 1)*dz   !z_trapezoid(2)
                 !print*,'
                 if(z1.lt.z_trap_temp+z_clearanceFactor*dz)then
                   i_createParticle = 0
                   i_trapezoid = n_trapezoid + 1
                 else
                   i_createParticle = 1
                 endif
               elseif(x1.ge.x_trapezoid(i_trapezoid,1))then
                 trapezoid_slope1 = 
     &           (z_trapezoid(i_trapezoid,2)-z_trapezoid(i_trapezoid,1))
     &          /(x_trapezoid(i_trapezoid,2)-x_trapezoid(i_trapezoid,1))            
                 z_trap_temp = z_trapezoid(i_trapezoid,1) + 
     &           (x1 - x_trapezoid(i_trapezoid,1))*trapezoid_slope1
                 !print*,'
                 if(z1.lt.z_trap_temp+z_clearanceFactor*dz)then
                   i_createParticle = 0
                   i_trapezoid = n_trapezoid + 1
                 else
                   i_createParticle = 1
                 endif
               else
                 i_createParticle = 1
               endif
              enddo   !End of: do while(i_createParticle.eq.0)
             endif
             !- Floating Bodies -  For now, All objects are assumed to be square cylinders horizontally orientated
             if(num_FB.gt.0.and.i_createParticle.eq.1)then
               i_num_FB = 0
               i_outsideObject = 1
               do while (i_num_FB.lt.num_FB.and.i_outsideObject.eq.1)
                 i_num_FB = i_num_FB + 1
                 !- Define clearance limits around object
                 if(iopt_FloatingObject_type(i_num_FB).eq.1)then   !Square cylinder
                   !- Transform particle coordinates to Object Frame of reference -
                   !- rotate about Centre of Gravity (CofG) - 
                   dx_newParticle = x1 - Box_XC(i_num_FB)
                   dz_newParticle = z1 - Box_ZC(i_num_FB)
                   !- transformed coordinates -
                   x1_Tr =   dx_newParticle*cos(body_Angle(i_num_FB))
     &                     + dz_newParticle*sin(body_Angle(i_num_FB))
                   z1_Tr = - dx_newParticle*sin(body_Angle(i_num_FB))
     &                     + dz_newParticle*cos(body_Angle(i_num_FB))
                   clearance_factor = 0.99
                   xtmin = - 0.5*XcylinderDimension(i_num_FB)
     &                     - clearance_factor*dx   !+ 0.5*dx
                   xtmax =   0.5*XcylinderDimension(i_num_FB)
     &                     + clearance_factor*dx   !+ 0.50*dx
                   ztmin = - 0.5*ZcylinderDimension(i_num_FB)
     &                     - clearance_factor*dz   !+ 0.5*dz
                   ztmax =   0.5*ZcylinderDimension(i_num_FB)
     &                     + clearance_factor*dz
                   if(x1_Tr.gt.xtmin.and.x1_Tr.lt.xtmax)then
                     if((z1_Tr.gt.ztmin.and.z1_Tr.lt.ztmax))then  !WP is within x-range of cylinder but above or below
                       i_outsideObject = 0
                     else
                       i_outsideObject = 1
                     endif
                   else
                     i_outsideObject = 1
                   end if
                 elseif(iopt_FloatingObject_type(i_num_FB).eq.2)then   !Antifer
                   !- Transform particle coordinates to Object Frame of reference -
                   !- rotate about bottom left-hand corner - 
                   dx_newParticle = x1 - x_cyl_min(i_num_FB)
                   dz_newParticle = z1 - z_cyl_min(i_num_FB)
                   !- transformed coordinates -
                   x1_Tr =   dx_newParticle*cos(body_Angle(i_num_FB))
     &                     + dz_newParticle*sin(body_Angle(i_num_FB))
                   z1_Tr = - dx_newParticle*sin(body_Angle(i_num_FB))
     &                     + dz_newParticle*cos(body_Angle(i_num_FB))
                   clearance_factor = 0.99
                   xtmin = !- 0.5*antiferBottom(i_num_FB)
     &                     - clearance_factor*dx   !+ 0.5*dx
                   xtmax =   antiferBottom(i_num_FB)
     &                     + clearance_factor*dx   !+ 0.50*dx
                   z_CoG_temp = Antifer_CoG_Z_ratio(i_num_FB)
     &                      *antiferHeight(i_num_FB)
                   ratio_temp = Antifer_HeightRatio(i_num_FB)
     &                         /Antifer_XRatio(i_num_FB)
                   b_temp = antiferBottom(i_num_FB)
                   if(x1_Tr.lt.0.5*b_temp)then
                     ztmaxo = min(
c     &                 (-z_CoG_temp + ratio_temp*(x1_Tr + 0.5*b_temp)),
c     &                 (antiferHeight(i_num_FB)-z_CoG_temp)   )
     &                  ratio_temp*(x1_Tr),
     &                 (antiferHeight(i_num_FB))   )
                     xtmin2 = x1_Tr - (z1_Tr - ztmaxo)/ratio_temp
c     &                     - clearance_factor*dx
                     xtmax2 = xtmax
                   else
                     ztmaxo = min(
c     &                 (-z_CoG_temp - ratio_temp*(x1_Tr - 0.5*b_temp)),
c     &                 (antiferHeight(i_num_FB)-z_CoG_temp)   )
     &                 ( - ratio_temp*(x1_Tr - b_temp)),
     &                 (antiferHeight(i_num_FB))   )
                     xtmin2 = xtmin
                     xtmax2 = x1_Tr + (z1_Tr - ztmaxo)/ratio_temp
c     &                      + clearance_factor*dx
                     if(ztmaxo.lt.0.0)then
                       ztmaxo = 0.0
                     endif
                   endif
                   ztmin =          - clearance_factor*dz
                   ztmax =   ztmaxo + clearance_factor*dz
                   !if(ztmax.lt.antiferHeight(i_num_FB))then
                   !  ztmax = ztmax + 2.0*clearance_factor*dz
                   !endif
                   if(x1_Tr.gt.xtmin.and.x1_Tr.lt.xtmax)then
                     if((z1_Tr.gt.ztmin .and.z1_Tr.lt.ztmax))then  !WP is within x-range of cylinder but above or below
                       i_outsideObject = 0
                     else
                       if(z1_Tr.lt.antiferHeight(i_num_FB)+0.5*dz)then
                         if(x1_Tr.gt.xtmin2.and.x1_Tr.lt.xtmax2)then
                           i_outsideObject = 0
                         else
                           i_outsideObject = 1
                         endif
                       else
                         i_outsideObject = 1
                       endif
                     endif
                   else
                     i_outsideObject = 1
                   end if
                 else
                   print*,'Invalid floating Object type '
                   print*,'in subroutine fill_part'
                   print*,'iopt_FloatingObject_type = ',
     &                     iopt_FloatingObject_type(i_num_FB)
                   stop
                 endif
               enddo
               if(i_outsideObject.eq.1)then
                 i_createParticle = 1
               else
                 i_createParticle = 0
               end if
             endif
                          
             if(i_createParticle.eq.1)then
             
c               if(iBC.eq.2)then   !Dalrymple BC
c                   nn=nn+1       
c                   call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
c                   call pressure(nn,ZZmax,dx,dy,dz,expont,g)
c               else if(iBC.eq.1)then   !Repulsive BC
                 hk2 = (x1-N_ini*dx)*tan(beta)+j*dy*tan(theta)   
                 if(z1.gt.(1.001*hk2+0.25*dz).and.
     &              z1.lt.1.001*ZZmax.and.
     &              x1.gt.0.999*XXmin_initial.and.
     &              x1.lt.0.999*XXmax_initial+0.5*dx)then
                   nn=nn+1       
                   call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
                   call pressure(nn,ZZmax,dx,dy,dz,expont,g)
                 end if
c               end if
               
               if (lattice.eq.2) then !BCC	
                
                   x1=XXmin+(i-N_ini+0.5)*dx
                   y1=0.
                   z1=ZZinf+(k-L_ini+0.5)*dz
                   hk2 = (x1-N_ini*dx)*tan(beta)+j*dy*tan(theta)
                   z1=ZZinf+(k-L_ini+0.5)*dz
                   if(x1.le.1.001*XXmax.and.z1.le.1.001*ZZmax) then ! 1.01 to avoid  
cc                     nn=nn+1                                    !roundoff
cc                     call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
cc                     call pressure(nn,ZZmax,dx,dy,dz,expont,g)
                     if(iBC.eq.2)then   !Dalrymple BC
                         nn=nn+1       
                         call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
                         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
                     else if(iBC.eq.1)then   !Repulsive BC
                       if(z1.gt.1.001*hk2+0.25*dz.and.
     &                     z1.lt.1.001*ZZmax.and.
     &                     x1.gt.0.999*XXmin_initial.and.
     &                     x1.lt.0.999*XXmax_initial+0.5*dx)then
                         nn=nn+1       
                         call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
                         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
                       endif
                     end if
                   endif
                endif
             end if   !End of:   if(i_createParticle.eq.1)then
           enddo
         enddo
       enddo
       
       return       
       end
       
c	______
      
c	______

c	____________________ SUBROUTINE DROP

      subroutine drop(nn,g,dx,dy,dz,expont)
      include 'common.gen2D'

      double precision Xcen,Zcen,Radius


	    write(*,*) ' Center of the drop: Xcen, Zcen ??'  
	    read (*,*)  Xcen, Zcen			
	    write(*,*)  Xcen, Zcen

	    write(*,*) ' Radius'
	    read (*,*)  Radius
	    write(*,*)  Radius

		write(*,*) ' Velocity in X, Z'   
	    read (*,*)  uveloc,wveloc
	    write(*,*)  uveloc,wveloc
          vveloc=0

	    XXmin=Xcen-Radius
	    XXmax=Xcen+Radius
	    YYmin=0
	    YYmax=0
          ZZmin=Zcen-Radius
	    ZZmax=Zcen+Radius

	    N_end=nint(XXmax/dx)
	    N_ini=nint(XXmin/dx)

	    M_end=0
	    M_ini=0

	    L_end=nint(ZZmax/dz)
	    ZZmax=L_end*dz
	    L_ini=nint(ZZmin/dz)
	    ZZmin=L_ini*dz
	    write(*,*) 'Recalculated max and min Z ',ZZmax,ZZmin

			    
	do i=N_ini,N_end
	do j=M_ini,M_end
	do k=L_ini,L_end
	       Xdif=(i+0.5)*dx-Xcen
	       Xdif=Xdif*Xdif
	       Ydif=0
           
	       Zdif=(k+0.5)*dz-Zcen
	       Zdif=Zdif*Zdif

             if(abs(sqrt(Xdif+Ydif+Zdif)).le.radius) then
                nn=nn+1
                call pos_veloc(nn,(i+0.5)*dx,(j+0.5)*dy,
     +			  (k+0.5)*dz,uveloc,vveloc,wveloc)
	          call pressure(nn,ZZmax,dx,dy,dz,expont,g)
	          if (lattice.eq.2) then !BCC
				  nn=nn+1
                 call pos_veloc(nn,(i+1)*dx,(j+1)*dy,(k+1)*dz,
     +	             uveloc,vveloc,wveloc)
                 call pressure(nn,ZZmax,dx,dy,dz,expont,g)
	          endif
             endif
	enddo
	enddo
	enddo

      return

	end

c	______

c	____________________ SUBROUTINE WALL

      subroutine wall(nn,dx,dy,dz,beta,theta)
      include 'common.gen2D'

      double precision Xwall,Zwall

	write(*,*) ' Wall position in X coordinates ??'
	write(*,*) ' X Wall ??'
	read(*,*) Xwall

	write(*,*) ' Wall height ??'
	write(*,*) ' Zwall ??'
	read(*,*) Zwall

	    YYmin=0.
	    YYmax=0.
	    M_ini=0
	    M_end=0
	
	L_end=nint(Zwall/dz)
c	L_ini=0

	write(*,*)M_end,M_ini,L_end,L_ini
	write(*,*)'Begining of wall at ',nn+1

	do j=M_ini,M_end
	    hk=Xwall*tan(beta)+j*dy*tan(theta)
	    L_ini=nint(hk/dz)
	    do k=L_ini,L_end
	         nn=nn+1
	     call pos_veloc(nn,Xwall,j*dy,k*dz,0.,0.,0.)
		 
	       if((j+0.5)*dy.le.YYmax.and.(k+0.5)*dz.le.Zwall) then
  	         nn=nn+1
	         call pos_veloc(nn,Xwall+0.5*dx,(j+0.5)*dy,
     +		                   (k+0.5)*dz,0.,0.,0.)
	       endif
	   enddo
	enddo

	write(*,*)'End of wall at ',nn

	return
	end


c      ______________________ SUBROUTINE OBSTACLE


      subroutine obstacle(nn,dx,dy,dz,beta,theta,iopt_obst)
      include 'common.gen2D'

      double precision x1,z1


      write(*,*) ' Which kind of obstacle'
      write(*,*) ' (1) Solid '
      write(*,*) ' (2) With Solid Walls'
      read(*,*) iopt_kind
      write(*,*) iopt_kind

      if(iBC.eq.1.and.iopt_kind.eq.1)then
        print*
      write(*,*) 'Invalid Option with Monaghan BC'
        stop
      endif
      
      write(*,*) ' Density of points'
      read(*,*) ndens
      write(*,*) ndens
      dx1=dx/ndens
      dy1=dy/ndens
      dz1=dz/ndens

      write(*,*) ' Cube containing particles'
      write(*,*) ' XMin, Xmax ??'
      read(*,*) XXmin,XXmax
      write(*,*) XXmin,XXmax

      YYmin=0.
      YYmax=0.

      
      write(*,*) ' ZMin, Zmax ??'
      read(*,*) ZZmin,ZZmax
      write(*,*) ZZmin,ZZmax

      write(*,*) ' slope in X direction??'
      read(*,*) slope
      write(*,*) slope
      slope=slope*pi/180

      if (slope.eq.pi/2) then 	
        valtan_inv=0
        Xtop=1.01*XXmax
      else 	
        valtan_inv =1./tan(slope)
        kmax=nint(ZZmax/dz1)
        Xtop=1.01*(XXmax+kmax*dx1*valtan_inv)
      endif

      Nstepsx=nint((XXmax-XXmin)/dx1)+1
      Nstepsy=1
      Nstepsz=nint((ZZmax-ZZmin)/dz1)+1

      write(*,*) 'Npnts x/y/z', nstepsx,nstepsy,nstepsz
      write(*,*) 'Beginning of obstacle at ',nn+1

      N_ini=nint(XXmin/dx1)
      N_end=nint(XXmax/dx1)
      i0=nint(vlx1/dx1)

      M_ini=0
      M_end=0

      if(iBC.eq.2)then


        do i=N_ini,N_end
        do j=1,Nstepsy

          hk=(i-i0)*dx1*tan(beta)+j*dy1*tan(theta)
          ZZinf=max(ZZmin,hk)

          L_ini=nint(ZZinf/dz1)
          L_end=nint(ZZmax/dz1)

          do k=L_ini,L_end
              if(iopt_kind.eq.1) then
                x1=XXmin+(i-N_ini)*dx1+(k-L_ini)*valtan_inv*dx1
                y1=YYmin+(j-1)*dy1
                z1=ZZinf+(k-L_ini)*dz1
                nn=nn+1
                call pos_veloc(nn,x1,y1,z1,0.,0.,0.)

                x1=XXmin+(i-N_ini-0.5)*dx1+(k-L_ini)*valtan_inv*dx1
                y1=YYmin+(j-0.5)*dy1
                z1=ZZinf+(k-L_ini-0.5)*dz1

                if(x1.le.Xtop.and.y1.le.1.01*YYmax.and.
     +            z1.le.1.01*ZZmax) then ! 1.01 to avoid the roundoff
                  nn=nn+1
                  call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
                endif
              else if(iopt_kind.eq.2) then
                if (i.eq.N_ini.or.i.eq.N_end.or.k.eq.L_ini
     +               .or.k.eq.L_end) then
                    x1=XXmin+(i-N_ini)*dx1+(k-L_ini)*valtan_inv*dx1
                    y1=YYmin+(j-1)*dy1
                    z1=ZZinf+(k-L_ini)*dz1
                    nn=nn+1
                    call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
   
                    x1=XXmin+(i-N_ini-0.5)*dx1+(k-L_ini)*valtan_inv*dx1
                    y1=YYmin+(j-0.5)*dy1
                    z1=ZZinf+(k-L_ini-0.5)*dz1
                    nn=nn+1
                    call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
                endif
              else
                write(*,*) ' Not a valid option'
                stop
              endif
            enddo
          enddo	
        enddo

      elseif(iBC.eq.1)then

c       -- Block face 1 -  x = const = N_ini*dx1
        i = N_ini
        do j = M_ini,M_end

          hk=(i-i0)*dx1*tan(beta)+j*dy1*tan(theta)
          ZZinf=max(ZZmin,hk)
          L_ini=nint(ZZinf/dz1)
          L_end=nint(ZZmax/dz1)

          do k = L_ini,L_end
            nn=nn+1
            x1=XXmin+(i-N_ini)*dx1+(k-L_ini)*valtan_inv*dx1
            y1 = j*dy1
            z1 = k*dz1
            call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
            if(k.eq.L_ini)then
              i_minus1 = nn
              i_plus1  = nn + 1
            else if(k.eq.L_end)then
              i_minus1 = nn - 1
              i_plus1  = nn
            else
              i_minus1 = nn - 1
              i_plus1  = nn + 1
            end if
             !-- Normal information - neighbour data for        --
             !-- Repulsive Boundary Particles (BPs)             --
             iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
             iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end do
        end do
        nn_final_Face1 = nn

c       -- Block face 2 -  x = const + width = N_end*dx1
        i = N_end
        do j = M_ini,M_end

          hk=(i-i0)*dx1*tan(beta)+j*dy1*tan(theta)
          ZZinf=max(ZZmin,hk)
          L_ini=nint(ZZinf/dz1)
          L_end=nint(ZZmax/dz1)

          do k = L_ini,L_end
            nn=nn+1
            x1=XXmin+(i-N_ini)*dx1+(k-L_ini)*valtan_inv*dx1
            y1 = j*dy1
            z1 = k*dz1
            call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
            if(k.eq.L_ini)then
              i_minus1 = nn + 1
              i_plus1  = nn
            else if(k.eq.L_end)then
              i_minus1 = nn
              i_plus1  = nn - 1
            else
              i_minus1 = nn + 1
              i_plus1  = nn - 1
            end if
             !-- Normal information - neighbour data for        --
             !-- Repulsive Boundary Particles (BPs)             --
             iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
             iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end do
        end do
        nn_final_Face2 = nn

c       -- Block face 3 - z = const = 
        L_end=nint(ZZmin/dz1)
        k = L_end
        do i = N_ini,N_end
          do j = M_ini,M_end
            nn=nn+1
            x1=XXmin+(i-N_ini)*dx1+(k-L_ini)*valtan_inv*dx1
            y1 = j*dy1
            z1 = k*dz1
            call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
            if(i.eq.N_ini)then
              i_minus1 = nn + 1
              i_plus1  = nn
            else if(i.eq.N_end)then
              i_minus1 = nn
              i_plus1  = nn - 1
            else
              i_minus1 = nn + 1
              i_plus1  = nn - 1
            end if
            !-- Normal information - neighbour data for        --
            !-- Repulsive Boundary Particles (BPs)             --
            iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
            iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end do
        end do
        nn_final_Face3 = nn

c       -- Block face 4 - z = const = 
        L_end=nint(ZZmax/dz1)
        k = L_end
        do i = N_ini,N_end
          do j = M_ini,M_end
            nn=nn+1
            x1=XXmin+(i-N_ini)*dx1+(k-L_ini)*valtan_inv*dx1
            y1 = j*dy1
            z1 = k*dz1
            call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
            if(i.eq.N_ini)then
              i_minus1 = nn
              i_plus1  = nn + 1
            else if(i.eq.N_end)then
              i_minus1 = nn - 1
              i_plus1  = nn
            else
              i_minus1 = nn - 1
              i_plus1  = nn + 1
            end if
            !-- Normal information - neighbour data for        --
            !-- Repulsive Boundary Particles (BPs)             --
            iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
            iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end do
        end do
        nn_final_Face4 = nn

      endif     !End of:   if(iBC.eq.2)then

      write(*,*) 'End of obstacle at ',nn

      write(55,*) iopt_obst
      write(55,*) XXmin
      write(55,*) XXmax
      write(55,*) YYmin
      write(55,*) YYmax
      write(55,*) ZZinf
      write(55,*) ZZmax
      write(55,*) slope

      return
      end

c      ______________________ SUBROUTINE TRAPEZOID


      subroutine trapezoid(nn,dx,dy,dz,beta,theta,iopt_obst)
      include 'common.gen2D'

      double precision x1,z1,x2,z2

       n_trapezoid = n_trapezoid + 1
       if(n_trapezoid_max.gt.n_trapezoid_max)then
         print*
         print*,'ERROR in SPHYSICSgen_2D.f'
         print*,'n_trapezoid.gt.n_trapezoid_max'
         print*,'Adjust n_trapezoid_max in common.gen2D'
         stop
       endif
       print*,'Trapezoid Number = ',n_trapezoid
       write(*,*) ' Enter (x,z)-start  of trapezoid'
       read(*,*) x_trapezoid(n_trapezoid,1),z_trapezoid(n_trapezoid,1)
       write(*,*) x_trapezoid(n_trapezoid,1),z_trapezoid(n_trapezoid,1)
       write(*,*) ' Enter (x,z)-start  of trapezoid top'
       read(*,*) x_trapezoid(n_trapezoid,2),z_trapezoid(n_trapezoid,2)
       write(*,*) x_trapezoid(n_trapezoid,2),z_trapezoid(n_trapezoid,2)
       write(*,*) ' Enter (x,z)-finish of trapezoid top'
       read(*,*) x_trapezoid(n_trapezoid,3),z_trapezoid(n_trapezoid,3)
       write(*,*) x_trapezoid(n_trapezoid,3),z_trapezoid(n_trapezoid,3)
       write(*,*) ' Enter (x,z)-finish of trapezoid'
       read(*,*) x_trapezoid(n_trapezoid,4),z_trapezoid(n_trapezoid,4)
       write(*,*) x_trapezoid(n_trapezoid,4),z_trapezoid(n_trapezoid,4)
       
c       -- Trapezoid face 1: Up-slope - 
        y1 = 0.0
        M_start  = 0
        M_finish = 0 
        N_start  = nint(x_trapezoid(n_trapezoid,1)/dx) !+ 1
        N_finish = nint(x_trapezoid(n_trapezoid,2)/dx) !+ 1 
        L_start  = nint(z_trapezoid(n_trapezoid,1)/dz) !+ 1
        L_finish = nint(z_trapezoid(n_trapezoid,2)/dz) !+ 1
        trapezoid_slope1 = 
     &     (z_trapezoid(n_trapezoid,2) - z_trapezoid(n_trapezoid,1))/
     &     (x_trapezoid(n_trapezoid,2) - x_trapezoid(n_trapezoid,1))   
        ntemp = nn
        do i = N_start,N_finish
          nn=nn+1
          x1 = i*dx
          z1 = z_trapezoid(n_trapezoid,1)
     &       + (x1 - x_trapezoid(n_trapezoid,1))*trapezoid_slope1
          call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
          if(iBC.eq.1)then
            if(i.eq.N_start)then
              i_minus1 = nn 
              i_plus1  = nn + 1
            else if(i.eq.N_finish)then
              i_minus1 = nn - 1 
              i_plus1  = nn
            else
              i_minus1 = nn - 1
              i_plus1  = nn + 1
            end if
            !-- Normal information - neighbour data for        --
            !-- Repulsive Boundary Particles (BPs)             --
            iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
            iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end if
          if(iBC.ne.1)then
            x2 = x1 + 0.5*dx
            y2 = y1 + 0.5*dy     
            z2 = z1 + 0.5*dz
            nn=nn+1
            call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
          end if
          if(i.eq.N_start)then
            print*,'Start of Trap, x1,z1 ',x1,z1
          endif
        end do
        z_lastPoint = z1
        if(iBC.eq.1)then
          !-- Normal information - neighbour data for              --
          !-- first corner particle on cylinder side               --
          iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
          iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
          !-- last corner particle on cylinder side                --
          iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
          iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
          ntemp = nn
        end if

c       -- Trapezoid face 2: Flat-Top - 
        y1 = 0.0
        M_start  = 0
        M_finish = 0 
        N_start  = nint(x_trapezoid(n_trapezoid,2)/dx) !+ 1
        N_finish = nint(x_trapezoid(n_trapezoid,3)/dx) !+ 1 
        z1 = z_lastPoint
        ntemp = nn
        do i = N_start,N_finish
          nn=nn+1
          x1 = i*dx
          call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
          if(iBC.eq.1)then
            if(i.eq.N_start)then
              i_minus1 = nn 
              i_plus1  = nn + 1
            else if(i.eq.N_finish)then
              i_minus1 = nn - 1 
              i_plus1  = nn
            else
              i_minus1 = nn - 1
              i_plus1  = nn + 1
            end if
            !-- Normal information - neighbour data for        --
            !-- Repulsive Boundary Particles (BPs)             --
            iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
            iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end if
          if(iBC.ne.1)then
            x2 = x1 + 0.5*dx
            y2 = y1 + 0.5*dy     
            z2 = z1 + 0.5*dz
            nn=nn+1
            call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
          end if
        end do
        x_TopFinish = x1
        if(iBC.eq.1)then
          !-- Normal information - neighbour data for              --
          !-- first corner particle on cylinder side               --
          iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
          iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
          !-- last corner particle on cylinder side                --
          iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
          iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
          ntemp = nn
        end if

c       -- Trapezoid face 3: Down-slope - 
        y1 = 0.0
        M_start  = 0
        M_finish = 0 
        N_start  = nint(x_trapezoid(n_trapezoid,3)/dx) !+ 1
        N_finish = nint(x_trapezoid(n_trapezoid,4)/dx) !+ 1 
        L_start  = nint(z_lastPoint) !+ 1   !int(z_trapezoid(3)/dz) + 1
        L_finish = nint(z_trapezoid(n_trapezoid,4)/dz) !+ 1
        trapezoid_slope1 =
     &     (z_trapezoid(n_trapezoid,4) - z_trapezoid(n_trapezoid,3))/
     &     (x_trapezoid(n_trapezoid,4) - x_trapezoid(n_trapezoid,3))   
        print*,'x_TopFinish ',x_TopFinish
        print*,'N_start*dx ',N_start*dx
        !if(N_start*dx.lt.x_TopFinish) N_start = N_start + 1
        ntemp = nn
        do i = N_start,N_finish
          nn=nn+1
          x1 = i*dx
          !z = z_trapezoid(3) + (x - x_trapezoid(3))*trapezoid_slope1
          z1 = z_lastPoint + (x1 - x_TopFinish)*trapezoid_slope1
          call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
          if(iBC.eq.1)then
            if(i.eq.N_start)then
              i_minus1 = nn 
              i_plus1  = nn + 1
            else if(i.eq.N_finish)then
              i_minus1 = nn - 1 
              i_plus1  = nn
            else
              i_minus1 = nn - 1
              i_plus1  = nn + 1
            end if
            !-- Normal information - neighbour data for        --
            !-- Repulsive Boundary Particles (BPs)             --
            iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
            iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          end if
          if(iBC.ne.1)then
            x2 = x1 + 0.5*dx
            y2 = y1 + 0.5*dy     
            z2 = z1 + 0.5*dz
            nn=nn+1
            call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
          end if
        end do
        if(iBC.eq.1)then
          !-- Normal information - neighbour data for              --
          !-- first corner particle on cylinder side               --
          iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
          iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
          !-- last corner particle on cylinder side                --
          iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
          iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
          ntemp = nn
        end if
c       -- End of Trapezoid
        print*,'End of Trap, x1,z1 ',x1,z1


       write(*,*) 'End of trapezoid, n_trapezoid= ',n_trapezoid,
     &            ' at nn=',nn

      if(n_trapezoid.eq.1)then

      endif

       return
       end


c	____________________ SUBROUTINE GATE

      subroutine gate(nn,dx,dy,dz,beta,theta)
      include 'common.gen2D'

      double precision Xwall

	write(*,*) ' Gate position in X coordinates ??'
	write(*,*) ' X Wall ??'
	read(*,*) Xwall
	write(*,*) Xwall
      
	YYmin=0.
	YYmax=0.

	write(*,*) ' Start and end of the gate in Z ??'
	write(*,*) ' ZZMin, ZZmax ??'
	read(*,*) ZZmin,ZZmax
	write(*,*) ZZmin,ZZmax

	write(*,*) ' GATE Velocity  ?? '
	write(*,*) ' Vx,Vz??'
	read(*,*) VXgate,VZgate
	write(*,*) VXgate,VZgate

	write(*,*) ' Initial time to start gate movement??'
	write(*,*) ' tgate ??'
	read(*,*) tgate
	write(*,*) tgate

      M_end=0
      M_ini=0

	L_end=nint(ZZmax/dz)
	L_ini=nint(ZZmin/dz)


	write(*,*) M_end,M_ini,L_end,L_ini
	write(*,*)'Begining of gate at ',nn+1

	ngate_ini=nn+1

	do j=M_ini,M_end
	   do k=L_ini,L_end
	         nn=nn+1
	      call pos_veloc(nn,Xwall,j*dy,k*dz,0.,0.,0.)

	       if((j+0.5)*dy.le.YYmax.and.(k+0.5)*dz.le.ZZmax) then
  	        nn=nn+1
	        call pos_veloc(nn,Xwall+0.5*dx,(j+0.5)*dy,
     +                                (k+0.5)*dz,0.,0.,0.)
	       endif
	   enddo
	enddo

	write(*,*)'End of gate at ',nn

	ngate_end=nn

	write(77,*) ngate_ini
	write(77,*) ngate_end
	write(77,*) VXgate,VZgate
      write(77,*) tgate

	return
	end


c	____________________ SUBROUTINE WAVEMAKER

      subroutine wavemaker(nn,dx,dy,dz,beta,theta)
      include 'common.gen2D'

      double precision XX_extra
       
       character (LEN=40) :: paddle_fileName 
       
       print*
       print*,'--- Generating Wavemaker ---'
       print*
       print*, ' Enter Paddle Type  ??'
       print*,' 1: Piston '
       print*,' 2: Piston-Flap (see Schafer et al. 1996)'
       print*,' 3: Piston with prescribed motion'
       read(*,*)  i_PaddleType
       print*,    i_PaddleType
       if(i_PaddleType.ne.1.and.i_PaddleType.ne.2
     &                     .and.i_PaddleType.ne.3)then
         print*,'Incorrect Paddle Type'
         stop
       end if
       
       print*,' Wavemaker Centre position in X coordinates ??'
       print*,' (Generally X_PaddleCentre = 0.5*stroke)'
       read(*,*)  X_PaddleCentre
       print*,    X_PaddleCentre
       
       if(i_paddleType.eq.2)then  !Info Needed for piston-flap paddle
         print*,' Info Needed for piston-flap paddle '
         print*,' Enter paddle Still Water Level (SWL)'
         read(*,*) paddle_SWL
         print*,   paddle_SWL
         print*,' Enter piston-flap flap_length'
         print*,' (distance under the bed of the pivot point' 
         read(*,*) flap_length
         print*,   flap_length
         paddle_fileName = 'noPaddleFileName'
       else if(i_paddleType.eq.3)then  !Info Needed for piston paddle with presribed motion
         print*
         print*,' Info Needed for piston-flap paddle '
         print*,' Enter filename of prescribed motion'
         print*,
     &        ' NOTE file format must be: time, x-position, u-velocity'
         read(*,*) paddle_fileName       
         print*,' paddle_fileName = ',paddle_fileName
         !- Read initial position from file -
         open(61,file=paddle_fileName)
           read(61,*)time_temp,X_PaddleStart_3,u_temp
         close(61)
         print*
         print*,'X_PaddleStart_3 = ',X_PaddleStart_3
         print*
         print*,' Enter paddle Still Water Level (SWL)'
         read(*,*) paddle_SWL
         print*,   paddle_SWL
         flap_length = 0.0
       else
         paddle_SWL  = 0.0
         flap_length = 0.0
         paddle_fileName = 'noPaddleFileName'
       end if
       
       YYmin=0.
       YYmax=0.
        
       write(*,*) ' Start and end of the wavemaker in Z ??'
       write(*,*) ' ZZMin, ZZmax ??'
       read(*,*)  ZZmin,ZZmax
       write(*,*) ZZmin,ZZmax
       
       write(*,*) ' Initial time of wavemaker ??'
       write(*,*) ' twavemaker ??'
       read(*,*)  twavemaker
       write(*,*) twavemaker
       
       write(*,*) ' Number of frequencies  ?? '
       read(*,*)  Nfreq
       write(*,*) Nfreq
       
     
       do n=1,Nfreq
       
         if(i_paddleType.eq.2)then
           write(*,*) ' Wavemaker Stroke = 2*Amplitude   ?? '
           read(*,*)  stroke
           write(*,*) stroke
         elseif(i_paddleType.eq.3)then
           write(*,*) ' Enter ZERO Wavemaker Stroke = 2*Amplitude   ?? '
           read(*,*)  stroke
           write(*,*) stroke
         elseif(i_paddleType.eq.1)then
           write(*,*) ' Amplitude   ?? '
           read(*,*)  A_wavemaker(n)
           write(*,*) A_wavemaker(n)
	     stroke=2*A_wavemaker(n)
         end if

         write(*,*) ' Period  ?? '
         read(*,*)  Period(n)
         write(*,*) Period(n)
     
         write(*,*) ' Phase  ?? '
         read(*,*)  Phase(n)
         write(*,*) Phase(n)
       
         write(*,*) ' twinitial  ?? ' 
         read(*,*)  twinitial(n)
         write(*,*) twinitial(n)

       enddo
                    
       M_start  = 0
       M_finish = 0
       
       L_end=nint(ZZmax/dz)
       L_ini=nint(ZZmin/dz)
       L_start = L_ini
       L_finish = L_end
       
       
       print*,'M_start, M_finish ',M_start, M_finish
       print*,'L_start, L_finish ',L_start, L_finish
       write(*,*)'Start of wavemaker at ',nn+1
       
       nwavemaker_ini=nn+1
       
       X_PaddleStart = 0.5*stroke
       
      if(iBC.eq.2)then        !Dalrymple BC - Fixed Fluid Particles

        do j=M_start,M_finish
          do k=L_start,L_finish
            if(k*dz.gt.ZZmin)then
              if(i_paddleType.eq.2)then
                XX_extra =  X_PaddleCentre 
     &          +0.5*(1.0+(k*dz - paddle_SWL)/
     &          (flap_length+paddle_SWL))*stroke
              elseif(i_paddleType.eq.3)then
                 XX_extra  = X_PaddleStart_3
              else
                 XX_extra  = X_PaddleCentre
              end if
              nn=nn+1
              call pos_veloc(nn,XX_extra,
     +                             j*dy,k*dz,0.,0.,0.)       
            end if
            if((j+0.5)*dy.le.YYmax.and.(k+0.5)*dz.le.ZZmax) then
              nn=nn+1
              call pos_veloc(nn,XX_extra+0.5*dx,
     &                           (j+0.5)*dy,(k+0.5)*dz,0.,0.,0.)
            endif
              
          enddo
        enddo
        
      else if(iBC.eq.1)then   !Repulsive BC - Repulsive Force Particles with Normals
        
       i=0
       do k=L_start,L_finish
         do j=M_start,M_finish
           !-- Start point is full extension --
           if(i_paddleType.eq.2)then
             XX_extra =  X_PaddleCentre 
     &       +0.5*(1.0+(k*dz - paddle_SWL)/
     &       (flap_length+paddle_SWL))*stroke
           elseif(i_paddleType.eq.3)then
             XX_extra  = X_PaddleStart_3
           else
             XX_extra  = X_PaddleCentre
           end if
           nn=nn+1
           call pos_veloc(nn,XX_extra,j*dy,k*dz,0.,0.,0.)
            !2-D Normals
            if(k.eq.L_start)then
               i_minus1 = nn + 1
               i_plus1  = nn 
             else if(k.eq.L_finish)then
               i_minus1 = nn 
               i_plus1  = nn - 1
             else
               i_minus1 = nn + 1
               i_plus1  = nn - 1
             end if
             !-- Normal information - neighbour data for        --
             !-- Repulsive Boundary Particles (BPs)             --
             iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
             iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour

         enddo
       enddo
      
      end if     !End of:  if(iBC.eq.1)then
       
      write(*,*)'  End of wavemaker at ',nn
       
       
       nwavemaker_end=nn
       
       Xwavemaker=X_PaddleCentre
       
       write(66,*) i_paddleType
       write(66,*) nwavemaker_ini
       write(66,*) nwavemaker_end
       write(66,*) X_PaddleCentre
       write(66,*) X_PaddleStart
       write(66,*) paddle_SWL
       write(66,*) flap_length
       write(66,*) stroke
       write(66,*) twavemaker
       write(66,*) paddle_fileName
       write(66,*) Nfreq
       
       do n=1,Nfreq       
         write(66,*) A_wavemaker(n)
         write(66,*) Period(n)
         write(66,*) Phase(n)
         write(66,*) twinitial(n)       
       enddo 
    
       close(66)
    
       return
       end


c     ____
c      ____________________ SUBROUTINE RaichlenWedge_Particles

      subroutine RaichlenWedge_Particles(nn,N1,N,M,dx,dy,dz,beta)
      include 'common.gen2D'

      double precision x1,z1,x2,z2
      double precision dxW,dyW,dzW

c     This test is based loosely on the experiments of Synolakis et al. (####)

c       -- Read in Wedge parameters --
        write(*,*) 'Enter block-top elevation above SWL'
        read(*,*) bigDelta
        write(*,*) bigDelta
        write(*,*) 'Enter specific Weight'
        read(*,*) specificWeight
        write(*,*) specificWeight
        write(*,*) 'Enter block_length, block_height, block_width'
        read(*,*) block_length, block_height, block_width
        write(*,*) block_length, block_height, block_width
        
        bslope = tan(beta)
        x_beach_start = (int((vlx1/dx))+1)*dx
        shoreline_position = x_beach_start + h_SWL/bslope
        block_top = h_SWL + bigDelta
        block_xstart_point  = shoreline_position - 
     &                        ((block_height - bigDelta)/bslope)
        block_xfinish_point = block_xstart_point + block_length
        block_ystart_point  = 0.5*M*dy - 0.5*block_width
        block_yfinish_point = 0.5*M*dy + 0.5*block_width
        block_zstart_point  = (block_xstart_point-x_beach_start)*bslope
        if(iBC.eq.1)then
          write(*,*) 'Enter block discretisation size'
          read(*,*) dxW, dzW
          write(*,*) dxW, dzW
        else        
          dxW = dx
          dyW = dy
          dzW = dz
        end if
        print*
        print*,'-- Generating Sliding Wedge --'
        print*,'N1 ',N1
        print*,'N1*dx ',N1*dx
        print*,'x_beach_start ',x_beach_start
        print*,'bslope ',bslope
        print*,'h_SWL  ',h_SWL
        print*,'shoreline_position',shoreline_position
        print*,'block_top',block_top
        print*,'block_xstart_point ',block_xstart_point
        print*,'block_xfinish_point',block_xfinish_point
        print*,'block_ystart_point' ,block_ystart_point
        print*,'block_yfinish_point',block_yfinish_point
        print*,'block_zstart_point' ,block_zstart_point
        print*,'dxW',dxW
        print*,'dzW',dzW
        print*

        nn_initial = nn
 
c       -- Block face 3 - x = const = block_xstart_point
        x1 = block_xstart_point 
        y1 = block_ystart_point
        z1 = block_zstart_point
        x1 = real(int(x1/dxW) + 1)*dxW
        M_start  = 0
        M_finish = 0 
        L_start  = int(block_zstart_point /dzW)+1
        L_finish = int(block_top/dzW)+1
        do j = M_start,M_finish
          do k = L_start,L_finish
            nn=nn+1
            y1 = j*dyW
            z1 = k*dzW
            call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
            if(iBC.eq.1)then
                if(k.eq.L_start)then
                  i_minus1 = nn 
                  i_plus1  = nn + 1
                else if(k.eq.L_finish)then
                  i_minus1 = nn - 1 
                  i_plus1  = nn
                else
                  i_minus1 = nn - 1
                  i_plus1  = nn + 1
                end if
                !-- Normal information - neighbour data for        --
                !-- Repulsive Boundary Particles (BPs)             --
                iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
                iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
            end if
            !call pos_veloc(nn,x,y,z,0.,0.,0.)
            if(iBC.ne.1)then
              x2 = x1 + 0.5*dxW
              y2 = y1 + 0.5*dyW     
              z2 = z1 + 0.5*dzW
              nn=nn+1
              call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
            end if
            !z = z + dzW
          end do
        end do

c       -- Block face 4 - z = const = SWL + bigDelta
        x1 = block_xstart_point
        y1 = block_ystart_point
        z1 = (int(block_top/dzW)+1)*dzW
        M_start  = 0
        M_finish = 0 
        N_start  = int(block_xstart_point /dxW)+1
        N_finish = int(block_xfinish_point/dxW)+1
        do i = N_start,N_finish
          do j = M_start,M_finish
            nn=nn+1
            x1 = i*dxW
            y1 = j*dyW
            call pos_veloc(nn,x1,y1,z1,0.,0.,0.)
            if(iBC.eq.1)then
                if(k.eq.N_start)then
                  i_minus1 = nn 
                  i_plus1  = nn + 1
                else if(k.eq.N_finish)then
                  i_minus1 = nn - 1 
                  i_plus1  = nn
                else
                  i_minus1 = nn - 1
                  i_plus1  = nn + 1
                end if
                !-- Normal information - neighbour data for        --
                !-- Repulsive Boundary Particles (BPs)             --
                iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
                iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
            end if
            if(iBC.ne.1)then
              x2 = x1 + 0.5*dxW
              y2 = y1 + 0.5*dyW     
              z2 = z1 + 0.5*dzW
              nn=nn+1
              call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
            end if
          end do
        end do
c       -- End of Wedge

       write(*,*) 'no. of Wedge particles: ',nn
      
      write(88,*)bslope
111   format(3e18.10)
      open(unit=87,
     &  file='Benchmark_4_initialPositions.dat',
     &     status='unknown')
      do i = nn_initial+1,nn
        write(87,111)xp(i),zp(i)
      end do
      close(87)

      return
      end

c      ____________________ SUBROUTINE FloatingBodies_Particles

      subroutine FloatingBodies_Particles(nn,N1,N,M,dx,dy,dz,beta)
      include 'common.gen2D'

      write(*,*)
      write(*,*) ' Which type of Floating Object' 
      write(*,*) ' (1) Square '
      write(*,*) ' (2) Antifer'
      write(*,*) ' (3) Tetrapod (In 2-D a Tripod)'
      read(*,*) iopt_FloatingObject_temp
      write(*,*) iopt_FloatingObject_temp

      if(iopt_FloatingObject_temp.eq.1)then
        call FloatingSquare_Particles(nn,N1,N,M,dx,dy,dz,beta)
      elseif(iopt_FloatingObject_temp.eq.2)then
        call Antifer_Particles(nn,N1,N,M,dx,dy,dz,beta)
      elseif(iopt_FloatingObject_temp.eq.4)then
        print*
        print*,'ERROR: Invalid Option'
        print*,'You will need to write your own subroutine '
        print*,'for a new shape'
        stop
        !call OtherShape_Particles(nn,N1,N,M,dx,dy,dz,beta)
      else
        print*
        write(*,*) 'Invalid Option '
        write(*,*) 'Not yet coded'
        stop
      endif

      end subroutine
c      ____________________ SUBROUTINE FloatingSquare_Particles

      subroutine FloatingSquare_Particles(nn,N1,N,M,dx,dy,dz,beta)
      include 'common.gen2D'

      double precision x1,z1,x2,z2


      num_FB = num_FB + 1
      if(num_FB.gt.num_FB_max)then
        print*,'Number of floating Bodies exceeds max value'
        print*,'num_FB.gt.num_FB_max'
        print*,'Adjust num_FB_max in common.gen2D, num_FB_max = ',
     &          num_FB_max
        stop
      endif
      iopt_FloatingObject_type(num_FB) = 1
c     -- Read in Floating Body parameters --
      write(*,*) 'Enter square cylinder X & Z Dimensions'
      read(*,*)  XcylinderDimension(num_FB),ZcylinderDimension(num_FB)      != 0.1*2.0
      write(*,*) XcylinderDimension(num_FB),ZcylinderDimension(num_FB)     
      print*
      write(*,*) 'Enter specific Weight'
      read(*,*)  FB_specificWeight(num_FB)    != 1.0
      write(*,*) FB_specificWeight(num_FB)
      print*
      write(*,*) 'Enter x,z of Bottom Left Corner'
      read(*,*)  x_BottomLeft, z_BottomLeft
      write(*,*) x_BottomLeft, z_BottomLeft
      print*
      write(*,*) 'Enter (x,z) shift of Centre of Gravity '
      write(*,*) 'for use in Parallel Axis Theorem '
      read(*,*)  xdiff_CofG, zdiff_CofG  !0.0*dz   !Shift Box Centre of Gravity to make slightly unbalanced
      write(*,*) xdiff_CofG, zdiff_CofG
      print*
      write(*,*) 'Enter initial U,W velocity of Object'
      read(*,*)  bigU(num_FB), bigW(num_FB)
      write(*,*) bigU(num_FB), bigW(num_FB)
      print*
      write(*,*)
     &     'Enter initial Body Angle (deg) and Rotation Rate (Omega)'
      read(*,*)  body_Angle(num_FB), bigOmega(num_FB)       !Positive anticlockwise
      write(*,*) body_Angle(num_FB), bigOmega(num_FB)
      !- Converting from deg to radians -
      if(body_Angle(num_FB).gt.180.0)then
        body_Angle(num_FB) = body_Angle(num_FB) - 360.0
      endif
      body_Angle(num_FB) = body_Angle(num_FB)*pi/180.0

      print*
      write(*,*) 'Enter coefficient of Friction'
      read(*,*)  friction_coeff(num_FB)
      write(*,*) friction_coeff(num_FB)

      print*

       cylinderDensity(num_FB) = rho0* FB_specificWeight(num_FB)
       bigMass(num_FB) = cylinderDensity(num_FB)
     &          *(XcylinderDimension(num_FB)*ZcylinderDimension(num_FB))
       bigInertiaYY(num_FB) = bigMass(num_FB)
     &          *(XcylinderDimension(num_FB)**2
     &          + ZcylinderDimension(num_FB)**2)/12.0
       Box_XC(num_FB) = x_BottomLeft + 0.5*XcylinderDimension(num_FB)
c     &                               - 0.5*dx + dx 
       Box_ZC(num_FB) = z_BottomLeft + 0.5*ZcylinderDimension(num_FB)
c     &                               + dz + 0.5*dz
       Box_XC(num_FB) = Box_XC(num_FB) + xdiff_CofG
       Box_ZC(num_FB) = Box_ZC(num_FB) + zdiff_CofG
       bigInertiaYY(num_FB) = bigInertiaYY(num_FB)
     &                + bigMass(num_FB)
     &                * (xdiff_CofG**2 + zdiff_CofG**2)  !Parallel Axis Theroem
       print*
       print*,'-- Generating Floating Body --'
       print*,'Floating Body Number : ',num_FB
       print*,'Initial Conditions:'
       print*,'XcylinderDimension   ',XcylinderDimension(num_FB)
       print*,'ZcylinderDimension   ',ZcylinderDimension(num_FB)
       print*,'cylinderDensity      ',cylinderDensity(num_FB)
       print*,'bigMass              ',bigMass(num_FB)
       print*,'bigInertiaYY         ',bigInertiaYY(num_FB)
       print*,'xdiff_CofG           ',xdiff_CofG
       print*,'zdiff_CofG           ',zdiff_CofG
       print*,'Box_XC, Box_ZC       ',Box_XC(num_FB), Box_ZC(num_FB)
       print*,'bigU, bigW           ',bigU(num_FB), bigW(num_FB)
       print*,'body_Angle, bigOmega ',
     &         body_Angle(num_FB), bigOmega(num_FB)

       nn_initial = nn

       x_cyl_min(num_FB) = x_BottomLeft
       x_cyl_max(num_FB) = x_BottomLeft + XcylinderDimension(num_FB)
       z_cyl_min(num_FB) = z_BottomLeft
       z_cyl_max(num_FB) = z_BottomLeft + ZcylinderDimension(num_FB)
       print*,'x_cyl_min, x_cyl_max    ',
     &         x_cyl_min(num_FB), x_cyl_max(num_FB)
       print*,'z_cyl_min, z_cyl_max    ',
     &         z_cyl_min(num_FB), z_cyl_max(num_FB)
       M_start  = 0
       M_finish = 0 
       !- Be careful here, sometimes nint does not give correct answer -
       !- int should be used                                           -
       N_start  = nint( x_cyl_min(num_FB)/dx)+1
       N_finish = nint((x_cyl_min(num_FB)
     &                 + XcylinderDimension(num_FB))/dx)+1
       L_start  = nint( z_cyl_min(num_FB)/dz)+1
       L_finish = nint((z_cyl_min(num_FB)
     &                 + ZcylinderDimension(num_FB))/dz)+1

       !print*,'x_cyl_min(num_FB)/dx        ',x_cyl_min(num_FB)/dx
       !print*,'int( x_cyl_min(num_FB)/dx)  ',int( x_cyl_min(num_FB)/dx)
       !print*,'nint( x_cyl_min(num_FB)/dx) ',nint( x_cyl_min(num_FB)/dx)
       print*,'N_start, N_finish       ',N_start, N_finish
       print*,'L_start, L_finish       ',L_start, L_finish

       !2-D, therefore no y-component
       y1 = 0.
       v1 = 0.
 
       !x = 0.5*(xb_max - xb_min) - 0.5*cylinderDimension - 130.0*dx
       !x = x + 10.0*dx  ! + 11.0*dx!Asymmetrical box position
       !y = -0.5*cylinderDimension - 2.0*dy + 0.5*dy !+ 12.0*dy
       
       x = (N_start - 1) * dx
       z = (L_start - 1) * dz

       if(z.lt.z_cyl_min(num_FB))then
         print*
         print*,'Adjusting z_min of object'
         z = z_cyl_min(num_FB)
         print*
       endif

       z = z - dz
       
       print*,'Initial x,z ',x,z
       !stop

       !-- Side 1:  x = x1
       !print*,'Side 1:  x = x1'
       ntemp = nn
       do k = L_start,L_finish
         z=z+dz
         nn=nn+1
         dx_Box = x - Box_XC(num_FB)
         dz_Box = z - Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=Box_XC(num_FB) + xp_temp
         z1=Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(k.eq.L_start)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(k.eq.L_finish)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
         !print*,'k, x1,z1 ',k,x1,z1
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
         ntemp = nn
       end if
       
       !-- Side 2:  z = z2
       !print*,'Side 2:  z = z2'
       x = x - dx
       do i = N_start,N_finish
         x=x+dx
         nn=nn+1
         dx_Box = x - Box_XC(num_FB)
         dz_Box = z - Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=Box_XC(num_FB) + xp_temp
         z1=Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(i.eq.N_start)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(i.eq.N_finish)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
         ntemp = nn
       end if
       
       !-- Side 3:  x = x2
       !print*,'Side 3:  x = x2'
       z = z + dz
       do k = L_start,L_finish
         z=z-dz
         nn=nn+1
         dx_Box = x - Box_XC(num_FB)
         dz_Box = z - Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=Box_XC(num_FB) + xp_temp
         z1=Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(k.eq.L_start)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(k.eq.L_finish)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
         ntemp = nn
       end if
       !-- Side 4:  z = z1
       !print*,'Side 4:  z = z1'
       x = x + dx
       do i = N_start,N_finish
         x=x-dx
         nn=nn+1
         dx_Box = x - Box_XC(num_FB)
         dz_Box = z - Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=Box_XC(num_FB) + xp_temp
         z1=Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(i.eq.N_start)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(i.eq.N_finish)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
       end if



       nb_FB(num_FB) = nn
       write(*,*) 'Object number: ',num_FB
       write(*,*) 'no. of Object particles: nb_FB',nb_FB(num_FB)

       write(88,*)num_FB
       write(88,*)bigMass(num_FB)
       write(88,*)bigInertiaYY(num_FB)
       write(88,*)XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
       write(88,*)cylinderDensity(num_FB)
       write(88,*)FB_SpecificWeight(num_FB)
       write(88,*)friction_coeff(num_FB)
       write(88,*)Box_XC(num_FB),Box_ZC(num_FB)
       write(88,*)bigU(num_FB),bigW(num_FB),bigOmega(num_FB)
       write(88,*)nb_FB(num_FB)

      return
      end

c      ____________________ SUBROUTINE Antifer_Particles

      subroutine Antifer_Particles(nn,N1,N,M,dx,dy,dz,beta)
      include 'common.gen2D'

      double precision x1,z1,x2,z2


      num_FB = num_FB + 1
      if(num_FB.gt.num_FB_max)then
        print*,'Number of floating Bodies exceeds max value'
        print*,'num_FB.gt.num_FB_max'
        print*,'Adjust num_FB_max in common.gen2D, num_FB_max = ',
     &          num_FB_max
        stop
      endif
      iopt_FloatingObject_type(num_FB) = 2
c     -- Read in Floating Body parameters --
      write(*,*) 'Enter square Antifer Dimensions'
      write(*,*) 'Lengths: Bottom, Top, Height'
      read(*,*)  antiferBottom(num_FB),
     &           antiferTop(num_FB),
     &           antiferHeight(num_FB)
      write(*,*) antiferBottom(num_FB),
     &           antiferTop(num_FB),
     &           antiferHeight(num_FB)
      !-- Calculate Antifer X & Z dimensions
      Antifer_HeightRatio(num_FB) = antiferHeight(num_FB)
     &                             /antiferBottom(num_FB)
      Antifer_XRatio(num_FB) = 0.5*
     &                      (1.0 - (antiferTop(num_FB))
     &                             /antiferBottom(num_FB))
      XcylinderDimension(num_FB) = antiferBottom(num_FB)
      ZcylinderDimension(num_FB) = antiferHeight(num_FB)
      write(*,*) XcylinderDimension(num_FB),ZcylinderDimension(num_FB)     
      print*
      write(*,*) 'Enter specific Weight'
      read(*,*)  FB_specificWeight(num_FB)    != 1.0
      write(*,*) FB_specificWeight(num_FB)
      print*
      write(*,*) 'Enter x,z of Bottom Left Corner'
      read(*,*)  x_BottomLeft, z_BottomLeft
      write(*,*) x_BottomLeft, z_BottomLeft
      print*
      write(*,*) 'Enter (x,z) shift of Centre of Gravity '
      write(*,*) 'for use in Parallel Axis Theorem '
      read(*,*)  xdiff_CofG, zdiff_CofG  !0.0*dz   !Shift Box Centre of Gravity to make slightly unbalanced
      write(*,*) xdiff_CofG, zdiff_CofG
      print*
      write(*,*) 'Enter initial U,W velocity of Object'
      read(*,*)  bigU(num_FB), bigW(num_FB)
      write(*,*) bigU(num_FB), bigW(num_FB)
      print*
      write(*,*)
     &     'Enter initial Body Angle (deg) and Rotation Rate (Omega)'
      read(*,*)  body_Angle(num_FB), bigOmega(num_FB)       !Positive anticlockwise
      write(*,*) body_Angle(num_FB), bigOmega(num_FB)
      !- Converting from deg to radians -
      if(body_Angle(num_FB).gt.180.0)then
        body_Angle(num_FB) = body_Angle(num_FB) - 360.0
      endif
      body_Angle(num_FB) = body_Angle(num_FB)*pi/180.0

      print*
      write(*,*) 'Enter coefficient of Friction'
      read(*,*)  friction_coeff(num_FB)
      write(*,*) friction_coeff(num_FB)

      print*

       x_cyl_min(num_FB) = x_BottomLeft
       x_cyl_max(num_FB) = x_BottomLeft + antiferBottom(num_FB)
       z_cyl_min(num_FB) = z_BottomLeft
       z_cyl_max(num_FB) = z_BottomLeft + antiferHeight(num_FB)
       print*,'x_cyl_min, x_cyl_max    ',
     &         x_cyl_min(num_FB), x_cyl_max(num_FB)
       print*,'z_cyl_min, z_cyl_max    ',
     &         z_cyl_min(num_FB), z_cyl_max(num_FB)
       M_start  = 0
       M_finish = 0 
       !- Be careful here, sometimes nint does not give correct answer -
       !- int should be used                                           -
       N_start_top  = nint((x_cyl_min(num_FB)
     &                    + Antifer_XRatio(num_FB)
     &                    *antiferBottom(num_FB))/dx)  !+1
       N_finish_top = nint((x_cyl_min(num_FB)
     &                    + XcylinderDimension(num_FB)
     &                  - Antifer_XRatio(num_FB)
     &                    *antiferBottom(num_FB) )/dx)  !+1
       N_start_bottom  = nint( x_cyl_min(num_FB)/dx)  !+1
       N_finish_bottom = nint((x_cyl_min(num_FB)
     &                 + antiferBottom(num_FB))/dx) !+1
       L_start = nint( z_cyl_min(num_FB)/dz)  !+1
       !if((z_cyl_min(num_FB)-L_start*dz).gt.0.5*dz)then
       !  L_start  = L_start + 1
       !endif
       L_finish = nint((z_cyl_min(num_FB)
     &                 + antiferHeight(num_FB))/dz) !+1

       !print*,'x_cyl_min(num_FB)/dx        ',x_cyl_min(num_FB)/dx
       !print*,'int( x_cyl_min(num_FB)/dx)  ',int( x_cyl_min(num_FB)/dx)
       !print*,'nint( x_cyl_min(num_FB)/dx) ',nint( x_cyl_min(num_FB)/dx)
       print*,'N_start_top, N_finish_top        ',
     &         N_start_top, N_finish_top
       print*,'N_start_bottom, N_finish_bottom  ',
     &         N_start_bottom, N_finish_bottom
       print*,'L_start, L_finish       ',
     &         L_start, L_finish

       !- Correcting antifer sizes to be integer dx & dz -
       print*
       print*,'Correcting antifer sizes to be integer dx & dz'
       antiferBottom(num_FB) = (N_finish_bottom - N_start_bottom)*dx
       antiferTop(num_FB) = (N_finish_top - N_start_top)*dx
       antiferHeight(num_FB) = (L_finish - L_start)*dz
       x_cyl_min(num_FB) = N_start_bottom*dx
       x_cyl_max(num_FB) = N_finish_bottom*dx
       z_cyl_min(num_FB) = L_start*dz
       z_cyl_max(num_FB) = L_finish*dz
       print*,'x_cyl_min, x_cyl_max    ',
     &         x_cyl_min(num_FB), x_cyl_max(num_FB)
       print*,'z_cyl_min, z_cyl_max    ',
     &         z_cyl_min(num_FB), z_cyl_max(num_FB)
       print*

       cylinderDensity(num_FB) = rho0* FB_specificWeight(num_FB)
       bigMass(num_FB) = cylinderDensity(num_FB)
     &   *0.5*antiferHeight(num_FB)
     &   *(antiferBottom(num_FB)+antiferTop(num_FB))
       a_temp = antiferTop(num_FB)
       b_temp = antiferBottom(num_FB)
       h_temp = antiferHeight(num_FB)
       a2_temp = a_temp**2
       b2_temp = b_temp**2
       h2_temp = h_temp**2
       ab_temp = a_temp*b_temp
       !For moment of inertia of symmetric trapezoid see 
       !http://www.efunda.com/math/areas/IsosTrapazoid.cfm
       bigInertiaYY(num_FB) = cylinderDensity(num_FB)
     &   *(16.*h2_temp*ab_temp + 4.*h2_temp*b2_temp
     &    + 4.*h2_temp*a2_temp + 3.*a2_temp*a2_temp
     &    + 6.*a2_temp*b2_temp + 6.*ab_temp*a_temp
     &    + 6.*ab_temp*b_temp  + 3.*b2_temp*b2_temp)
     &    /(144.0*(a_temp+b_temp))
c       bigInertiaYY(num_FB) = bigMass(num_FB)
c     &          *(XcylinderDimension(num_FB)**2)/12.0
       Antifer_CoG_Z_ratio(num_FB) = (2.0*antiferTop(num_FB)
     &              + antiferBottom(num_FB))
     &          /(3.0*(antiferBottom(num_FB)
     &                +antiferTop(num_FB)))
       x_CoG_temp = 0.5*antiferBottom(num_FB)
       z_CoG_temp = antiferHeight(num_FB)*Antifer_CoG_Z_ratio(num_FB)
c       Box_XC(num_FB) = x_BottomLeft+0.5*antiferBottom(num_FB)
cc     &                               - 0.5*dx + dx 
c       Box_ZC(num_FB) = z_BottomLeft
c     &    + antiferHeight(num_FB)*Antifer_CoG_Z_ratio(num_FB)
       Box_XC(num_FB) = x_cyl_min(num_FB)
     &                + x_CoG_temp*cos(body_Angle(num_FB))
     &                - z_CoG_temp*sin(body_Angle(num_FB))
c     &                               - 0.5*dx + dx 
       Box_ZC(num_FB) = z_cyl_min(num_FB)
     &                + x_CoG_temp*sin(body_Angle(num_FB))
     &                + z_CoG_temp*cos(body_Angle(num_FB))
c     &                               + dz + 0.5*dz
       Box_XC(num_FB) = Box_XC(num_FB) + xdiff_CofG
       Box_ZC(num_FB) = Box_ZC(num_FB) + zdiff_CofG
       bigInertiaYY(num_FB) = bigInertiaYY(num_FB)
     &                + bigMass(num_FB)
     &                * (xdiff_CofG**2 + zdiff_CofG**2)  !Parallel Axis Theroem
       print*
       print*,'-- Generating Antifer --'
       print*,'Floating Body Number : ',num_FB
       print*,'Initial Conditions:'
       print*,'antiferBottom  ',antiferBottom(num_FB)
       print*,'antiferTop     ',antiferTop(num_FB)
       print*,'antiferHeight  ',antiferHeight(num_FB)
       print*,'cylinderDensity         ',cylinderDensity(num_FB)
       print*,'bigMass                 ',bigMass(num_FB)
       print*,'bigInertiaYY            ',bigInertiaYY(num_FB)
       print*,'xdiff_CofG              ',xdiff_CofG
       print*,'zdiff_CofG              ',zdiff_CofG
       print*,'Box_XC, Box_ZC          ',Box_XC(num_FB), Box_ZC(num_FB)
       print*,'bigU, bigW              ',bigU(num_FB), bigW(num_FB)
       print*,'body_Angle, bigOmega    ',
     &         body_Angle(num_FB), bigOmega(num_FB)
       print*,'x_CoG_temp, z_CoG_temp ',x_CoG_temp, z_CoG_temp
       
       nn_initial = nn

       !2-D, therefore no y-component
       y1 = 0.
       v1 = 0.
 
       !x = 0.5*(xb_max - xb_min) - 0.5*cylinderDimension - 130.0*dx
       !x = x + 10.0*dx  ! + 11.0*dx!Asymmetrical box position
       !y = -0.5*cylinderDimension - 2.0*dy + 0.5*dy !+ 12.0*dy
       
       x = (N_start_bottom ) * dx
     &     - Antifer_XRatio(num_FB)*antiferBottom(num_FB)
     &     /(L_finish - L_start)
       z = (L_start - 1) * dz

       !if(z.lt.z_cyl_min(num_FB))then
       !  print*
       !  print*,'Adjusting z_min of object'
       !  z = z_cyl_min(num_FB)
       !  print*
       !endif

       !z = z - dz
       
       print*,'Initial x,z ',x,z
       !stop

       !-- Side 1:  x = x1
       !print*,'Side 1:  x = x1'
       ntemp = nn
       do k = L_start,L_finish
         x=x+Antifer_XRatio(num_FB)*antiferBottom(num_FB)
     &     /(L_finish - L_start)
         !z = k*dz   
         z=z+dz
         nn=nn+1
         dx_Box = x - x_cyl_min(num_FB)  !Box_XC(num_FB)
         dz_Box = z - z_cyl_min(num_FB)  !Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=x_cyl_min(num_FB) + xp_temp   !Box_XC(num_FB) + xp_temp
         z1=z_cyl_min(num_FB) + zp_temp    !Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(k.eq.L_start)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(k.eq.L_finish)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
         !print*,'k, x1,z1 ',k,x1,z1
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
         ntemp = nn
       end if
       
       !-- Side 2:  z = z2
       !print*,'Side 2:  z = z2'
       x = x - dx
       do i = N_start_top,N_finish_top
         !x = i*dx   
         x=x+dx
c     &     *(antiferTop(num_FB)/antiferBottom(num_FB))
         nn=nn+1
         dx_Box = x - x_cyl_min(num_FB)  !Box_XC(num_FB)
         dz_Box = z - z_cyl_min(num_FB)  !Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=x_cyl_min(num_FB) + xp_temp   !Box_XC(num_FB) + xp_temp
         z1=z_cyl_min(num_FB) + zp_temp    !Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(i.eq.N_start_top)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(i.eq.N_finish_top)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
         ntemp = nn
       end if
       
       !-- Side 3:  x = x2
       !print*,'Side 3:  x = x2'
       x=x-Antifer_XRatio(num_FB)*antiferBottom(num_FB)
     &     /(L_finish - L_start)
       z = z + dz
       do k = L_finish,L_start,-1
         x=x+Antifer_XRatio(num_FB)*antiferBottom(num_FB)
     &     /(L_finish - L_start)
         !z = k*dz   
         z=z-dz
         nn=nn+1
         dx_Box = x - x_cyl_min(num_FB)  !Box_XC(num_FB)
         dz_Box = z - z_cyl_min(num_FB)  !Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=x_cyl_min(num_FB) + xp_temp   !Box_XC(num_FB) + xp_temp
         z1=z_cyl_min(num_FB) + zp_temp    !Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(k.eq.L_start)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(k.eq.L_finish)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour     
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour     
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
         ntemp = nn
       end if
       !-- Side 4:  z = z1
       !print*,'Side 4:  z = z1'
       x = x + dx
       do i = N_finish_bottom,N_start_bottom,-1
        !x = i*dx   
        x=x-dx
        if(x.gt.x_cyl_min(num_FB))then
         nn=nn+1
         dx_Box = x - x_cyl_min(num_FB)  !Box_XC(num_FB)
         dz_Box = z - z_cyl_min(num_FB)  !Box_ZC(num_FB)
         xp_temp = dx_Box*cos(body_Angle(num_FB))
     &           - dz_Box*sin(body_Angle(num_FB))
         zp_temp = dx_Box*sin(body_Angle(num_FB))
     &           + dz_Box*cos(body_Angle(num_FB))
         x1=x_cyl_min(num_FB) + xp_temp   !Box_XC(num_FB) + xp_temp
         z1=z_cyl_min(num_FB) + zp_temp    !Box_ZC(num_FB) + zp_temp
         u1=bigU(num_FB) - bigOmega(num_FB)*(z1 - Box_ZC(num_FB))
         w1=bigW(num_FB) + bigOmega(num_FB)*(x1 - Box_XC(num_FB))
         call pos_veloc(nn,x1,y1,z1,u1,v1,w1)
         call pressure(nn,ZZmax,dx,dy,dz,expont,g)
         if(iBC.eq.1)then
           if(i.eq.N_start_bottom)then
             i_minus1 = nn 
             i_plus1  = nn + 1
           else if(i.eq.N_finish_bottom)then
             i_minus1 = nn - 1 
             i_plus1  = nn
           else
             i_minus1 = nn - 1
             i_plus1  = nn + 1
           end if
           !-- Normal information - neighbour data for        --
           !-- Repulsive Boundary Particles (BPs)             --
           iBP_Pointer_Info(nn,3) = i_minus1       !i-1 neighbour     
           iBP_Pointer_Info(nn,4) = i_plus1        !i+1 neighbour
          endif   !End of:   if(x.gt.x_cyl_min(num_FB))then
         end if
c         if(iBC.ne.1)then
c           x2 = x1 + 0.5*dx
c           y2 = y1 + 0.5*dy
c           z2 = z1 + 0.5*dz
c           nn=nn+1
c           call pos_veloc(nn,x2,y2,z2,0.,0.,0.)
c         end if
       end do
       if(iBC.eq.1)then
         !-- Normal information - neighbour data for              --
         !-- first corner particle on cylinder side               --
         iBP_Pointer_Info(ntemp+1,3) = ntemp+1        !i-1 neighbour
         iBP_Pointer_Info(ntemp+1,4) = ntemp+2        !i+1 neighbour
         !-- last corner particle on cylinder side                --
         iBP_Pointer_Info(nn,3) = nn -1         !i-1 neighbour
         iBP_Pointer_Info(nn,4) = nn            !i+1 neighbour
       end if

       !if(num_FB.eq.3)stop


       nb_FB(num_FB) = nn
       write(*,*) 'Object number: ',num_FB
       write(*,*) 'no. of Object particles: nb_FB',nb_FB(num_FB)

       write(88,*)num_FB
       write(88,*)bigMass(num_FB)
       write(88,*)bigInertiaYY(num_FB)
       write(88,*)XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
       write(88,*)cylinderDensity(num_FB)
       write(88,*)FB_SpecificWeight(num_FB)
       write(88,*)friction_coeff(num_FB)
       write(88,*)Box_XC(num_FB),Box_ZC(num_FB)
       write(88,*)bigU(num_FB),bigW(num_FB),bigOmega(num_FB)
       write(88,*)nb_FB(num_FB)

      return
      end

c      ____________________ SUBROUTINE Tetrapod_Particles

      subroutine Tetrapod_Particles(nn,N1,N,M,dx,dy,dz,beta)
      include 'common.gen2D'

        print*
        write(*,*) 'Invalid Option '
        write(*,*) ' Not yet coded'
        stop

      end subroutine
c     ____________________ SUBROUTINE WAVE

      subroutine wave(nn,XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax,
     +	    		     g,dx,dy,dz,expont)
      include 'common.gen2D'

      double precision x1,z1,hh

      write(*,*) ' Cube containing particles'
      write(*,*) ' XMin, Xmax ??'
      read(*,*) XXmin,XXmax
      write(*,*) XXmin,XXmax

      YYMin=0.
	YYmax=0.  

      read(*,*) YYmin,YYmax
      write(*,*) YYmin,YYmax
      write(*,*) ' ZZMin ?? '
      read(*,*) ZZmin
      write(*,*) ZZmin
      
      write(*,*) ' depth ??'
      read(*,*) D
      write(*,*) D
      write(*,*) ' amplitude ??'
      read(*,*) A
      write(*,*) A

      vlong=sqrt(0.75*(D**3)/A)

      write(*,*) ' x0 ??'
      read(*,*) x0
      write(*,*) x0

      Nstepsy=nint((YYmax-YYmin)/dy)+1
	Nstepsy=1

      N_end=nint(XXmax/dx)
      N_ini=nint(XXmin/dx)


      do j=1,Nstepsy	
        do i=N_ini,N_end
       
          x=i*dx
          alt=A/(cosh((x-x0)/vlong)*cosh((x-x0)/vlong))
          vx=alt*(sqrt(9.81/D))
          hh=D+alt
          L_ini=nint(ZZmin/dz)
          L_end=nint(h/dz)
          
          do K=L_ini,L_end
            x1=XXmin+(i-N_ini)*dx
            y1=YYmin+(j-1)*dy
            y1=0
            z1=ZZmin+(k-L_ini)*dz
            nn=nn+1
            call pos_veloc(nn,x1,y1,z1,vx,0.,0.)
            call pressure(nn,hh,dx,dy,dz,expont,g)
              

              x1=XXmin+(i-N_ini+0.5)*dx
              y1=0
              z1=ZZmin+(k-L_ini+0.5)*dz
              if(x1.le.1.01*XXmax.and.z1.le.1.01*ZZmax) then ! 1.01 to avoid the roundoff
                nn=nn+1
                call pos_veloc(nn,x1,y1,z1,vx,0.,0.)
                call pressure(nn,hh,dx,dy,dz,expont,g)
              endif

          enddo
        enddo
      enddo

      return
      end

c     ____

c     ____________________ SUBROUTINE NORMALS_CALC_2D
      
      subroutine normals_Calc_2D

      include 'common.gen2D'
      
      double precision dxpb1, dxpb2
      
      do i=1,nb
        
        xp_i = xp(i)
        zp_i = zp(i)
        if(iBP_Pointer_Info(i,3).eq.0.and.
     &     iBP_Pointer_Info(i,4).eq.0)then
          if(i.eq.1)then
            i_plus1  = i + 1
            i_minus1 = i
          else if(i.eq.nbf)then   !!Correct This For The Interface Between Two Processes!!
            i_plus1  = i
            i_minus1 = i - 1
          else if(i.eq.nbf+1)then
            i_plus1  = i + 1
            i_minus1 = i
          else if(i.eq.nb)then   !!Correct This For The Interface Between Two Processes!!
            i_plus1  = i
            i_minus1 = i - 1
          else
            i_plus1  = i + 1
            i_minus1 = i - 1
          end if
        else   !Prescribed special neighbours for difficult points
          i_minus1 = iBP_Pointer_Info(i,3)
          i_plus1  = iBP_Pointer_Info(i,4)
        end if
        
        if(i.eq.i_minus1)then
          xp_previous = xp(i) - (xp(i_plus1)-xp(i))
          zp_previous = zp(i) - (zp(i_plus1)-zp(i))
          xp_next = xp(i_plus1)
          zp_next = zp(i_plus1)
        else if(i.eq.i_plus1)then
          xp_previous = xp(i_minus1)
          zp_previous = zp(i_minus1)
          xp_next = xp(i) + (xp(i)-xp(i_minus1))
          zp_next = zp(i) + (zp(i)-zp(i_minus1))
        else
          xp_previous = xp(i_minus1)
          zp_previous = zp(i_minus1)
          xp_next = xp(i_plus1)
          zp_next = zp(i_plus1)
        end if
                
        iBP_Pointer_Info(i,3) = i_minus1
        iBP_Pointer_Info(i,4) = i_plus1
        BP_xz_Data(i,1) = xp(i)
        BP_xz_Data(i,2) = zp(i)
        
        dxn_forward=(xp_next-xp_i)
        dzn_forward=(zp_next-zp_i)
        dxn_back=(xp_i-xp_previous)
        dzn_back=(zp_i-zp_previous)
                
        rbn_back = sqrt(dxn_back*dxn_back + dzn_back*dzn_back)
        rbn_forward = sqrt(dxn_forward*dxn_forward +
     &                     dzn_forward*dzn_forward)
        dzn = dzn_forward + dzn_back
        dxn = dxn_forward + dxn_back
        rb=dxn*dxn+dzn*dzn
        xnb(i)=-dzn/sqrt(rb)
        znb(i)=dxn/sqrt(rb)

      enddo
      
      print*
      print*,'NORMALS calculated'
      print*

      end subroutine  
c     ____

c     ____________________ SUBROUTINE NORMALS_FILEWRITE_2D
c

      subroutine normals_FileWrite_2D
      
      include 'common.gen2D'

      
    
      open(21,file='NORMALS.init')

22    format(2(e16.8,' '),4(i6,' '),4(e16.8,' '))
      do i=1,nb
                
        
        iBP_Pointer_Info(i,1) = i  !Default MPI Rank
        iBP_Pointer_Info(i,2) = 0  !Default MPI Rank
        
        !-- Diagnostic Check --
        if(i.eq.0)then
          print*,'iBP_Pointer_Info(i,1)',
     &            iBP_Pointer_Info(i,1)     
          print*,'iBP_Pointer_Info(i,2)',
     &            iBP_Pointer_Info(i,2)     
          write(*,22) xnb(i),znb(i),
     &    iBP_Pointer_Info(i,1),
     &    iBP_Pointer_Info(i,2),
     &    iBP_Pointer_Info(i,3),
     &    iBP_Pointer_Info(i,4),
     &    BP_xz_Data(i,1),
     &    BP_xz_Data(i,2)
        end if
        !---------------------
        
        write(21,22) xnb(i),znb(i),
     &    iBP_Pointer_Info(i,1),
     &    iBP_Pointer_Info(i,2),
     &    iBP_Pointer_Info(i,3),
     &    iBP_Pointer_Info(i,4),
     &    BP_xz_Data(i,1),
     &    BP_xz_Data(i,2)
        
      enddo
      close(21)
      
      print*
      print*,'NORMALS.init file written'
      print*
      
      end subroutine  
c     ____

c     ____________________ SUBROUTINE PERIODICITY_CHECK
c

      subroutine periodicityCheck(drxp, dryp, drzp,
     &                            xp_i, xp_other,
     &                            yp_i, yp_other,
     &                            zp_i, zp_other)
      
      include 'common.gen2D'

      double precision dxpb1,dxpb2
      double precision dypb1,dypb2
      double precision dzpb1,dzpb2
      
      
c     ---  Subroutine modifies differences in position if periodicity is activated ---      
      
      
        !--- X-Direction Periodicity ---
        if(i_periodicOBs(1).eq.1)then
          if(drxp.gt.6.0*h)then
            dxpb1 = dble(xb_max+dpx)-dble(xp_other)   !xb_max_double-dble(xp_next)
            dxpb2 = dble(xp_i)-dble(xb_min)      !dble(xp_i)-xb_min_double
            drxp = -real(dxpb1+dxpb2)
          else if(drxp.lt.-6.0*h)then
            dxpb1 = dble(xb_max+dpx)-dble(xp_i)    !xb_max_double-dble(xp_i)
            dxpb2 = dble(xp_other)-dble(xb_min)     !dble(xp_next)-xb_min_double
            drxp = real(dxpb1+dxpb2)
          end if
        end if
        !--- Y-Direction Periodicity ---
        if(i_periodicOBs(2).eq.1)then
          if(dryp.gt.6.0*h)then
            dypb1 = dble(yb_max+dpy)-dble(yp_other)   !xb_max_double-dble(xp_next)
            dypb2 = dble(yp_i)-dble(yb_min)      !dble(xp_i)-xb_min_double
            dryp = -real(dypb1+dypb2)
          else if(dryp.lt.-6.0*h)then
            dypb1 = dble(yb_max+dpy)-dble(yp_i)    !xb_max_double-dble(xp_i)
            dypb2 = dble(yp_other)-dble(yb_min)     !dble(xp_next)-xb_min_double
            dryp = real(dypb1+dypb2)
          end if
        end if
        !--- Z-Direction Periodicity ---
        if(i_periodicOBs(3).eq.1)then
          if(drzp.gt.6.0*h)then
            dzpb1 = dble(zb_max+dpz)-dble(zp_other)   !xb_max_double-dble(xp_next)
            dzpb2 = dble(zp_i)-dble(zb_min)      !dble(xp_i)-xb_min_double
            drzp = -real(dzpb1+dzpb2)
          else if(drzp.lt.-6.0*h)then
            dzpb1 = dble(zb_max+dpz)-dble(zp_i)    !xb_max_double-dble(xp_i)
            dzpb2 = dble(zp_other)-dble(zb_min)     !dble(xp_next)-xb_min_double
            drzp = real(dzpb1+dzpb2)
          end if
        end if
      
      end subroutine  


c     ____________________ SUBROUTINE TOCOMPILE_gfortran

      subroutine tocompile_gfortran
      include 'common.gen2D'
      CHARACTER(LEN=10) :: FMT,FMT1
      CHARACTER(LEN=1)  :: TAB,LS
      
      TAB=CHAR(9)
      LS =CHAR(92)
     
      FMT="(A)"
      FMT1="(2A)"
      open(22,file='SPHYSICS.mak')
      write(22,FMT) 'FC=gfortran'
      write(22,FMT) 'OPTIONS= -O3'
      !- Uncomment as required -
c     &   //' -ftree-vectorize'
c     &   //' -ffast-math -funroll-loops'
c     &   //' -fbounds-check'
c     &   //' -frange-check'
c     &   //' -Wunderflow'
c     &   //' -Wuninitialized'
c     &   //' -ffpe-trap=invalid,zero,overflow'
c     &   //' -fstack-check'
c     &   //' -fstack-protector'
c     &   //' -ftrapv'
c     &   //' -ftracer'
c     &   //' -g'
c     &   //' -fimplicit-none'
c     &   //' -fno-automatic'
      write(22,FMT) 'srcdir=.'
      write(22,FMT) 'idir=../../execs'
      write(22,FMT) 'bakdir=../../execs.bak'
      write(22,FMT) 'objects=energy_2D.o recover_list_2D.o ini_divide_2D
     &.o \'
      write(22,FMT1)TAB,'keep_list_2D.o SPHYSICS_2D.o  \'
      write(22,FMT1)TAB,'getdata_2D.o check_limits_2D.o \'                  
      write(22,FMT1)TAB,'divide_2D.o vorticity_2D.o\'
      write(22,FMT1)TAB,'movingObjects_2D.o updateNormals_2D.o \'
      write(22,FMT1)TAB,'movingGate_2D.o \'
      write(22,FMT1)TAB,'movingPaddle_2D.o movingWedge_2D.o \'
      write(22,FMT1)TAB,'periodicityCorrection_2D.o \'

      !- Kernel Corrections
      if (i_kernelcorrection.eq.0) then
        write(22,FMT1)TAB,'ac_NONE_2D.o \'
        write(22,FMT1)TAB,'kernel_correction_NC_2D.o \'
      elseif (i_kernelcorrection.eq.2) then
        write(22,FMT1)TAB,'ac_KGC_2D.o \'
        write(22,FMT1)TAB,'kernel_correction_KGC_2D.o \'
        write(22,FMT1)TAB,'pre_self_KGC_2D.o \'
        write(22,FMT1)TAB,'pre_celij_KGC_2D.o \'
      elseif (i_kernelcorrection.eq.1) then
        write(22,FMT1)TAB,'ac_KC_2D.o \'
        write(22,FMT1)TAB,'kernel_correction_KC_2D.o \'
      endif

      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'ac_2D.o \'
        write(22,FMT1)TAB,'poute_2D.o \'
        write(22,FMT1)TAB,'gradients_calc_basic_2D.o \'
      else
        write(22,FMT1)TAB,'ac_Conservative_2D.o \'
        write(22,FMT1)TAB,'poute_Conservative_2D.o \'
        write(22,FMT1)TAB,'gradients_calc_Conservative_2D.o \'
      endif

      !- Boundary Conditions
      if(iBC.eq.1) then
        if(num_FB.eq.0)then
          write(22,FMT1)TAB,'monaghanBC_2D.o \'
c          write(22,FMT1)TAB,'monaghanBC_AngMomCorrection_2D.o \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,'self_BC_Monaghan_2D.o \'
            write(22,FMT1)TAB,'celij_BC_Monaghan_2D.o \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannConservative_2D.o \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannNonConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannNonConservative_2D.o \'
            endif
          endif
        else
          !- FLOATING BODIES
          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_2D.o \'
c          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_'
c     &                             //'AngMomCorrection_2D.o \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,
     &                'self_BC_Monaghan_FloatingBodies_2D.o \'
            write(22,FMT1)TAB,
     &               'celij_BC_Monaghan_FloatingBodies_2D.o \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.o \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.o \'
            endif
          endif
        endif
      else !(iBC.eq.2) 
        if(num_FB.eq.0)then
           write(22,FMT1)TAB,'self_BC_Dalrymple_2D.o \'
           write(22,FMT1)TAB,'celij_BC_Dalrymple_2D.o \'
        else
           write(22,FMT1)TAB,'self_BC_Dalrymple_FloatingBodies_2D.o \'
           write(22,FMT1)TAB,'celij_BC_Dalrymple_FloatingBodies_2D.o \'
        endif
      endif

      !- Rigid Body Motion 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'rigid_body_motion_2D.o \'
      else
        write(22,FMT1)TAB,'rigid_body_motion_Conservative_2D.o \'
      endif

      !- Variable Time Step 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'variable_time_step_2D.o \'
      else
        write(22,FMT1)TAB,'variable_time_step_Conservative_2D.o \'
      endif

      !- Approximate Riemann Solver
      if(iRiemannSolver.gt.0)then 
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'NonConservative_2D.o \'
        else
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'Conservative_2D.o \'
        endif
        i_limiterChoice = 1  !Temporary Default Value
        if(i_limiterChoice.eq.1)then
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.o \'
        else
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.o \'
        endif
      endif

      !- Viscosity
      if (i_viscos.eq.1) then
          write(22,FMT1)TAB,'viscosity_artificial_2D.o \'
      else if (i_viscos.eq.2) then
          write(22,FMT1)TAB,'viscosity_laminar_2D.o \'
      else
          write(22,FMT1)TAB,'viscosity_laminar+SPS_2D.o \'
      endif
		
      if(i_viscos.eq.3) then
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_SPS_2D.o \'
        else
          print*
          print*,'Option not yet available: '
          print*,'SPS model in Conservative formulation'
          stop
          write(22,FMT1)TAB,'correct_SPS_Conservative_2D.o \'
        endif
      else
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_2D.o \'
        else
          write(22,FMT1)TAB,'correct_Conservative_2D.o \'
        endif
      endif
        

c	    Kernels
c	    (1) Gaussian; (2)Quadratic;  (3)Cubic  ; (5) Wendland Quintic
      if (i_kernel.eq.1) then
             write(22,FMT1)TAB,'kernel_gaussian_2D.o \'
      else if (i_kernel.eq.2) then 
             write(22,FMT1)TAB,'kernel_quadratic_2D.o \'
      else if (i_kernel.eq.3) then 
             write(22,FMT1)TAB,'kernel_cubic_2D.o \'
      else if (i_kernel.eq.5) then
	     write(22,FMT1)TAB,'kernel_wendland5_2D.o \'
      endif

c         Equation of State
      if (i_EoS.eq.1) then
             write(22,FMT1)TAB,'EoS_Tait_2D.o \'
      else if (i_EoS.eq.2) then
             write(22,FMT1)TAB,'EoS_IdealGas_2D.o \'
      else if (i_EoS.eq.3) then
             write(22,FMT1)TAB,'EoS_Morris_2D.o \'
      else
             write(22,FMT1)TAB,'EoS_Poisson_2D.o \'
      endif

c         Density Reinitialisation
      if (i_densityFilter.eq.1) then
         write(22,FMT1)TAB,'densityFilter_Shepard_2D.o \'
         write(22,FMT1)TAB,'ac_Shepard_2D.o \'
         write(22,FMT1)TAB,'pre_celij_Shepard_2D.o \'
         write(22,FMT1)TAB,'pre_self_Shepard_2D.o \'
      else if (i_densityFilter.eq.2) then
         write(22,FMT1)TAB,'densityFilter_MLS_2D.o \'
         write(22,FMT1)TAB,'ac_MLS_2D.o \'
         write(22,FMT1)TAB,'LU_decomposition_2D.o \'
         write(22,FMT1)TAB,'pre_celij_MLS_2D.o \'
         write(22,FMT1)TAB,'pre_self_MLS_2D.o \'
      else   
         write(22,FMT1)TAB,'densityFilter_NONE_2D.o \'
      endif

c         Algorithm 
      if(i_ConservativeFormulation.eq.0)then
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_2D.o '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.o '
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_2D.o '
        else
          write(22,FMT1)TAB,'step_beeman_2D.o '
        endif
      elseif(i_ConservativeFormulation.eq.1)then   !ConservativeFormulation
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_'
     &                    //'Conservative_2D.o '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.o '
          print*
          print*,'Not a valid Option : '
          print*,'step_verlet in Conservative formulation'
          stop
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_Conservative_2D.o '
        else
          write(22,FMT1)TAB,'step_beeman_2D.o '
          print*
          print*,'Not a valid Option : '
          print*,'step_beeman in Conservative formulation'
          stop
        endif
      else
        print*
        print*,'Not a valid Option : '
        print*,'SPS model in Conservative formulation'
        stop
      endif


      write(22,FMT)'#'      
      write(22,FMT)'SPHYSICS_2D: $(objects)'
      write(22,FMT1)TAB,'$(FC) $(OPTIONS) -o SPHYSICS_2D $(objects)'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'if [ -d $(bakdir) ]; then \'
      write(22,FMT1)TAB,'echo "execs.bak Directory Exists"; else \'
      write(22,FMT1)TAB,'mkdir $(bakdir); \'
      write(22,FMT1)TAB,'echo "execs.bak Directory Created"; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'if [ -d $(idir) ]; then \'
      write(22,FMT1)TAB,'echo "execs Directory Exists"; else \'
      write(22,FMT1)TAB,'mkdir $(idir); \'
      write(22,FMT1)TAB,'echo "execs Directory Created"; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#' 
      write(22,FMT1)TAB,'-if [ -f $(idir)/SPHYSICS_2D ]; then \'
      write(22,FMT1)TAB,'mv -f $(idir)/SPHYSICS_2D $(bakdir)/; \'
      write(22,FMT1)TAB,'echo Old SPHYSICS_2D moved to execs.bak'
     &//' from execs; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'mv SPHYSICS_2D $(idir)'
      write(22,FMT1)TAB,'echo New SPHYSICS_2D moved to execs'
      write(22,FMT)'#'
      write(22,FMT)'clean:'
      write(22,FMT1)TAB,'rm *.o'
      write(22,FMT1)TAB,'rm *~' 
      write(22,FMT)'#'
      write(22,FMT)'%.o: %.f'
      write(22,FMT1)TAB,'$(FC) $(OPTIONS) -c -o $@ $<'
       
      close(22)
      
      return
      end	
c     ____  

c     ____________________ SUBROUTINE TOCOMPILE_ifort

      subroutine tocompile_ifort
      include 'common.gen2D'
      CHARACTER(LEN=10) :: FMT,FMT1
      CHARACTER(LEN=1)  :: TAB,LS
      
      TAB=CHAR(9)
      LS =CHAR(92)
     
      FMT="(A)"
      FMT1="(2A)"
      open(22,file='SPHYSICS.mak')
      write(22,FMT) 'FC=ifort'
      write(22,FMT) 'OPTIONS= -O3'
      !- Uncomment as required -
     &   //' -ipo'
c     &   //' -g -debug all -check all'
c     &   //' -implicitnone'
c     &   //' -warn unused'
c     &   //' -fp-stack-check -heap-arrays'
c     &   //' -ftrapuv'
c     &   //' -check pointers'
c     &   //' -check bounds'
c     &   //' -check uninit'
c     &   //' -fpstkchk'
c     &   //' -traceback'
      write(22,FMT) 'srcdir=.'
      write(22,FMT) 'idir=../../execs'
      write(22,FMT) 'bakdir=../../execs.bak'
      write(22,FMT) 'objects=energy_2D.o recover_list_2D.o ini_divide_2D
     &.o \'
      write(22,FMT1)TAB,'keep_list_2D.o SPHYSICS_2D.o  \'
      write(22,FMT1)TAB,'getdata_2D.o check_limits_2D.o \'                  
      write(22,FMT1)TAB,'divide_2D.o vorticity_2D.o\'
      write(22,FMT1)TAB,'movingGate_2D.o \'
      write(22,FMT1)TAB,'movingObjects_2D.o updateNormals_2D.o \'
      write(22,FMT1)TAB,'movingPaddle_2D.o movingWedge_2D.o \'
      write(22,FMT1)TAB,'periodicityCorrection_2D.o \'

      !- Kernel Corrections
      if (i_kernelcorrection.eq.0) then
        write(22,FMT1)TAB,'ac_NONE_2D.o \'
        write(22,FMT1)TAB,'kernel_correction_NC_2D.o \'
      elseif (i_kernelcorrection.eq.1) then
        write(22,FMT1)TAB,'ac_KGC_2D.o \'
        write(22,FMT1)TAB,'kernel_correction_KGC_2D.o \'
        write(22,FMT1)TAB,'pre_self_KGC_2D.o \'
        write(22,FMT1)TAB,'pre_celij_KGC_2D.o \'
      elseif (i_kernelcorrection.eq.2) then
        write(22,FMT1)TAB,'ac_KC_2D.o \'
        write(22,FMT1)TAB,'kernel_correction_KC_2D.o \'
      endif

      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'ac_2D.o \'
        write(22,FMT1)TAB,'poute_2D.o \'
        write(22,FMT1)TAB,'gradients_calc_basic_2D.o \'
      else
        write(22,FMT1)TAB,'ac_Conservative_2D.o \'
        write(22,FMT1)TAB,'poute_Conservative_2D.o \'
        write(22,FMT1)TAB,'gradients_calc_Conservative_2D.o \'
      endif


      !- Boundary Conditions
      if(iBC.eq.1) then
        if(num_FB.eq.0)then
          write(22,FMT1)TAB,'monaghanBC_2D.o \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,'self_BC_Monaghan_2D.o \'
            write(22,FMT1)TAB,'celij_BC_Monaghan_2D.o \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannConservative_2D.o \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannNonConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannNonConservative_2D.o \'
            endif
          endif
        else
          !- FLOATING BODIES
          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_2D.o \'
c          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_'
c     &                             //'AngMomCorrection_2D.o \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,
     &                'self_BC_Monaghan_FloatingBodies_2D.o \'
            write(22,FMT1)TAB,
     &               'celij_BC_Monaghan_FloatingBodies_2D.o \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.o \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.o \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.o \'
            endif
          endif
        endif
      else !(iBC.eq.2) 
        if(num_FB.eq.0)then
           write(22,FMT1)TAB,'self_BC_Dalrymple_2D.o \'
           write(22,FMT1)TAB,'celij_BC_Dalrymple_2D.o \'
        else
           write(22,FMT1)TAB,'self_BC_Dalrymple_FloatingBodies_2D.o \'
           write(22,FMT1)TAB,'celij_BC_Dalrymple_FloatingBodies_2D.o \'
        endif
      endif

      !- Rigid Body Motion 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'rigid_body_motion_2D.o \'
      else
        write(22,FMT1)TAB,'rigid_body_motion_Conservative_2D.o \'
      endif

      !- Variable Time Step 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'variable_time_step_2D.o \'
      else
        write(22,FMT1)TAB,'variable_time_step_Conservative_2D.o \'
      endif

      !- Approximate Riemann Solver
      if(iRiemannSolver.gt.0)then 
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'NonConservative_2D.o \'
        else
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'Conservative_2D.o \'
        endif
        i_limiterChoice = 1  !Temporary Default Value
        if(i_limiterChoice.eq.1)then
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.o \'
        else
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.o \'
        endif
      endif

      !- Viscosity
      if (i_viscos.eq.1) then
          write(22,FMT1)TAB,'viscosity_artificial_2D.o \'
      else if (i_viscos.eq.2) then
          write(22,FMT1)TAB,'viscosity_laminar_2D.o \'
      else
          write(22,FMT1)TAB,'viscosity_laminar+SPS_2D.o \'
      endif
		
      if(i_viscos.eq.3) then
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_SPS_2D.o \'
        else
          print*
          print*,'Option not yet available: '
          print*,'SPS model in Conservative formulation'
          stop
          write(22,FMT1)TAB,'correct_SPS_Conservative_2D.o \'
        endif
      else
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_2D.o \'
        else
          write(22,FMT1)TAB,'correct_Conservative_2D.o \'
        endif
      endif
        

c	    Kernels
c	    (1) Gaussian; (2)Quadratic;  (3)Cubic  ; (5) Wendland Quintic
      if (i_kernel.eq.1) then
             write(22,FMT1)TAB,'kernel_gaussian_2D.o \'
      else if (i_kernel.eq.2) then 
             write(22,FMT1)TAB,'kernel_quadratic_2D.o \'
      else if (i_kernel.eq.3) then 
             write(22,FMT1)TAB,'kernel_cubic_2D.o \'
      else if (i_kernel.eq.5) then
	     write(22,FMT1)TAB,'kernel_wendland5_2D.o \'
      endif

c         Equation of State
      if (i_EoS.eq.1) then
             write(22,FMT1)TAB,'EoS_Tait_2D.o \'
      else if (i_EoS.eq.2) then
             write(22,FMT1)TAB,'EoS_IdealGas_2D.o \'
      else if (i_EoS.eq.3) then
             write(22,FMT1)TAB,'EoS_Morris_2D.o \'
      else
             write(22,FMT1)TAB,'EoS_Poisson_2D.o \'
      endif

c         Density Reinitialisation
      if (i_densityFilter.eq.1) then
         write(22,FMT1)TAB,'densityFilter_Shepard_2D.o \'
         write(22,FMT1)TAB,'ac_Shepard_2D.o \'
         write(22,FMT1)TAB,'pre_celij_Shepard_2D.o \'
         write(22,FMT1)TAB,'pre_self_Shepard_2D.o \'
      else if (i_densityFilter.eq.2) then
         write(22,FMT1)TAB,'densityFilter_MLS_2D.o \'
         write(22,FMT1)TAB,'ac_MLS_2D.o \'
         write(22,FMT1)TAB,'LU_decomposition_2D.o \'
         write(22,FMT1)TAB,'pre_celij_MLS_2D.o \'
         write(22,FMT1)TAB,'pre_self_MLS_2D.o \'
      else   
         write(22,FMT1)TAB,'densityFilter_NONE_2D.o \'
      endif

c         Algorithm 
      if(i_ConservativeFormulation.eq.0)then
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_2D.o '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.o '
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_2D.o '
        else
          write(22,FMT1)TAB,'step_beeman_2D.o '
        endif
      elseif(i_ConservativeFormulation.eq.1)then   !ConservativeFormulation
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_'
     &                    //'Conservative_2D.o '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.o '
          print*
          print*,'Not a valid Option : '
          print*,'step_verlet in Conservative formulation'
          stop
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_Conservative_2D.o '
        else
          write(22,FMT1)TAB,'step_beeman_2D.o '
          print*
          print*,'Not a valid Option : '
          print*,'step_beeman in Conservative formulation'
          stop
        endif
      else
        print*
        print*,'Not a valid Option : '
        print*,'SPS model in Conservative formulation'
        stop
      endif


      write(22,FMT)'#'      
      write(22,FMT)'SPHYSICS_2D: $(objects)'
      write(22,FMT1)TAB,'$(FC) $(OPTIONS) -o SPHYSICS_2D $(objects)'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'if [ -d $(bakdir) ]; then \'
      write(22,FMT1)TAB,'echo "execs.bak Directory Exists"; else \'
      write(22,FMT1)TAB,'mkdir $(bakdir); \'
      write(22,FMT1)TAB,'echo "execs.bak Directory Created"; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'if [ -d $(idir) ]; then \'
      write(22,FMT1)TAB,'echo "execs Directory Exists"; else \'
      write(22,FMT1)TAB,'mkdir $(idir); \'
      write(22,FMT1)TAB,'echo "execs Directory Created"; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#' 
      write(22,FMT1)TAB,'-if [ -f $(idir)/SPHYSICS_2D ]; then \'
      write(22,FMT1)TAB,'mv -f $(idir)/SPHYSICS_2D $(bakdir)/; \'
      write(22,FMT1)TAB,'echo Old SPHYSICS_2D moved to execs.bak'
     &//' from execs; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'mv SPHYSICS_2D $(idir)'
      write(22,FMT1)TAB,'echo New SPHYSICS_2D moved to execs'
      write(22,FMT)'#'
      write(22,FMT)'clean:'
      write(22,FMT1)TAB,'rm *.o'
      write(22,FMT1)TAB,'rm *~' 
      write(22,FMT)'#'
      write(22,FMT)'%.o: %.f'
      write(22,FMT1)TAB,'$(FC) $(OPTIONS) -c -o $@ $<'
       
      close(22)
      
      return
      end	
c     ____
  
       
c     ____________________ SUBROUTINE TOCOMPILE_WINDOWS


      subroutine tocompile_ftn95
      include 'common.gen2D'
      CHARACTER(LEN=10) :: FMT,FMT1
      CHARACTER(LEN=1)  :: TAB
      
      TAB=CHAR(9)
     
       print*,'Using SilverFrost FTN95 '
	 
      FMT="(A)"
      FMT1="(2A)"
      open(22,file='SPHYSICS.mak')
      write(22,FMT) 'OPTIONS= /OPTIMISE'
      !- Debugging options, Uncomment as required -
c     &   //' /CHECK'
      write(22,FMT)
      write(22,FMT) 'OBJFILES=energy_2D.obj recover_list_2D.obj \'
      write(22,FMT1)TAB,'ini_divide_2D.obj keep_list_2D.obj \'
      write(22,FMT1)TAB,'SPHYSICS_2D.obj getdata_2D.obj \'
      write(22,FMT1)TAB,'check_limits_2D.obj \'                    
      write(22,FMT1)TAB,'divide_2D.obj \'
      write(22,FMT1)TAB,'movingObjects_2D.obj movingGate_2D.obj\'
      write(22,FMT1)TAB,'movingPaddle_2D.obj movingWedge_2D.obj \'
      write(22,FMT1)TAB,'updateNormals_2D.obj vorticity_2D.obj\'
      write(22,FMT1)TAB,'periodicityCorrection_2D.obj \'


      !- Kernel Corrections
      if (i_kernelcorrection.eq.0) then
        write(22,FMT1)TAB,'ac_NONE_2D.obj \'
        write(22,FMT1)TAB,'kernel_correction_NC_2D.obj \'
      elseif (i_kernelcorrection.eq.2) then
        write(22,FMT1)TAB,'ac_KGC_2D.obj \'
        write(22,FMT1)TAB,'kernel_correction_KGC_2D.obj \'
        write(22,FMT1)TAB,'pre_self_KGC_2D.obj \'
        write(22,FMT1)TAB,'pre_celij_KGC_2D.obj \'
      elseif (i_kernelcorrection.eq.1) then
        write(22,FMT1)TAB,'ac_KC_2D.obj \'
        write(22,FMT1)TAB,'kernel_correction_KC_2D.obj \'
      endif

      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'ac_2D.obj \'
        write(22,FMT1)TAB,'poute_2D.obj \'
        write(22,FMT1)TAB,'gradients_calc_basic_2D.obj \'
      else
        write(22,FMT1)TAB,'ac_Conservative_2D.obj \'
        write(22,FMT1)TAB,'poute_Conservative_2D.obj \'
        write(22,FMT1)TAB,'gradients_calc_Conservative_2D.obj \'
      endif

      !- Boundary Conditions
      if(iBC.eq.1) then
        if(num_FB.eq.0)then
          write(22,FMT1)TAB,'monaghanBC_2D.obj \'
c          write(22,FMT1)TAB,'monaghanBC_AngMomCorrection_2D.obj \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,'self_BC_Monaghan_2D.obj \'
            write(22,FMT1)TAB,'celij_BC_Monaghan_2D.obj \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannConservative_2D.obj \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannNonConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannNonConservative_2D.obj \'
            endif
          endif
        else
          !- FLOATING BODIES
          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_2D.obj \'
c          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_'
c     &                             //'AngMomCorrection_2D.obj \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,
     &                'self_BC_Monaghan_FloatingBodies_2D.obj \'
            write(22,FMT1)TAB,
     &               'celij_BC_Monaghan_FloatingBodies_2D.obj \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.obj \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.obj \'
            endif
          endif
        endif
      else !(iBC.eq.2) 
        if(num_FB.eq.0)then
         write(22,FMT1)TAB,
     &        'self_BC_Dalrymple_2D.obj \'
         write(22,FMT1)TAB,
     &        'celij_BC_Dalrymple_2D.obj \'
        else
         write(22,FMT1)TAB,
     &        'self_BC_Dalrymple_FloatingBodies_2D.obj \'
         write(22,FMT1)TAB,
     &        'celij_BC_Dalrymple_FloatingBodies_2D.obj \'
        endif
      endif

      !- Rigid Body Motion 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'rigid_body_motion_2D.obj \'
      else
        write(22,FMT1)TAB,'rigid_body_motion_Conservative_2D.obj \'
      endif

      !- Variable Time Step 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'variable_time_step_2D.obj \'
      else
        write(22,FMT1)TAB,'variable_time_step_Conservative_2D.obj \'
      endif

      !- Approximate Riemann Solver
      if(iRiemannSolver.gt.0)then 
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'NonConservative_2D.obj \'
        else
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'Conservative_2D.obj \'
        endif
        i_limiterChoice = 1  !Temporary Default Value
        if(i_limiterChoice.eq.1)then
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.obj \'
        else
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.obj \'
        endif
      endif

      !- Viscosity
      if (i_viscos.eq.1) then
          write(22,FMT1)TAB,'viscosity_artificial_2D.obj \'
      else if (i_viscos.eq.2) then
          write(22,FMT1)TAB,'viscosity_laminar_2D.obj \'
      else
          write(22,FMT1)TAB,'viscosity_laminar+SPS_2D.obj \'
      endif
		
      if(i_viscos.eq.3) then
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_SPS_2D.obj \'
        else
          print*
          print*,'Option not yet available: '
          print*,'SPS model in Conservative formulation'
          stop
          write(22,FMT1)TAB,'correct_SPS_Conservative_2D.obj \'
        endif
      else
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_2D.obj \'
        else
          write(22,FMT1)TAB,'correct_Conservative_2D.obj \'
        endif
      endif
        

c	    Kernels
c	    (1) Gaussian; (2)Quadratic;  (3)Cubic  ; (5) Wendland Quintic
      if (i_kernel.eq.1) then
             write(22,FMT1)TAB,'kernel_gaussian_2D.obj \'
      else if (i_kernel.eq.2) then 
             write(22,FMT1)TAB,'kernel_quadratic_2D.obj \'
      else if (i_kernel.eq.3) then 
             write(22,FMT1)TAB,'kernel_cubic_2D.obj \'
      else if (i_kernel.eq.5) then
	     write(22,FMT1)TAB,'kernel_wendland5_2D.obj \'
      endif

c         Equation of State
      if (i_EoS.eq.1) then
             write(22,FMT1)TAB,'EoS_Tait_2D.obj \'
      else if (i_EoS.eq.2) then
             write(22,FMT1)TAB,'EoS_IdealGas_2D.obj \'
      else if (i_EoS.eq.3) then
             write(22,FMT1)TAB,'EoS_Morris_2D.obj \'
      else
             write(22,FMT1)TAB,'EoS_Poisson_2D.obj \'
      endif

c         Density Reinitialisation
      if (i_densityFilter.eq.1) then
         write(22,FMT1)TAB,'densityFilter_Shepard_2D.obj \'
         write(22,FMT1)TAB,'ac_Shepard_2D.obj \'
         write(22,FMT1)TAB,'pre_celij_Shepard_2D.obj \'
         write(22,FMT1)TAB,'pre_self_Shepard_2D.obj \'
      else if (i_densityFilter.eq.2) then
         write(22,FMT1)TAB,'densityFilter_MLS_2D.obj \'
         write(22,FMT1)TAB,'ac_MLS_2D.obj \'
         write(22,FMT1)TAB,'LU_decomposition_2D.obj \'
         write(22,FMT1)TAB,'pre_celij_MLS_2D.obj \'
         write(22,FMT1)TAB,'pre_self_MLS_2D.obj \'
      else   
         write(22,FMT1)TAB,'densityFilter_NONE_2D.obj \'
      endif

c         Algorithm 
      if(i_ConservativeFormulation.eq.0)then
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_2D.obj '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.obj '
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_2D.obj '
        else
          write(22,FMT1)TAB,'step_beeman_2D.obj '
        endif
      elseif(i_ConservativeFormulation.eq.1)then   !ConservativeFormulation
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_'
     &                    //'Conservative_2D.obj '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.obj '
          print*
          print*,'Not a valid Option : '
          print*,'step_verlet in Conservative formulation'
          stop
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_Conservative_2D.obj '
        else
          write(22,FMT1)TAB,'step_beeman_2D.obj '
          print*
          print*,'Not a valid Option : '
          print*,'step_beeman in Conservative formulation'
          stop
        endif
      else
        print*
        print*,'Not a valid Option : '
        print*,'SPS model in Conservative formulation'
        stop
      endif

      write(22,FMT)      
      write(22,FMT)'SPHYSICS_2D.exe: $(OBJFILES)'

      write(22,FMT)      
      write(22,FMT)'clean:'
      write(22,FMT1)TAB,'del *.obj'
       
      close(22)

      return
      end
	  
c     ____________________ SUBROUTINE TOCOMPILE_WIN_IFORT

      subroutine tocompile_win_ifort
      include 'common.gen2D'
      CHARACTER(LEN=10) :: FMT,FMT1
      CHARACTER(LEN=1)  :: TAB

      TAB=CHAR(9)

      print*,'Using Intel Fortran for Windows'

      FMT="(A)"
      FMT1="(2A)"
      open(22,file='SPHYSICS.mak')
      write(22,FMT) 'OPTIONS= /NOLOGO'
      write(22,FMT) 'COPTIONS= /03'

      write(22,FMT)
      write(22,FMT) 'OBJFILES=energy_2D.obj recover_list_2D.obj \'
      write(22,FMT1)TAB,'ini_divide_2D.obj keep_list_2D.obj \'      
      write(22,FMT1)TAB,'SPHYSICS_2D.obj getdata_2D.obj \'
      write(22,FMT1)TAB,'check_limits_2D.obj \'
      write(22,FMT1)TAB,'divide_2D.obj \'
      write(22,FMT1)TAB,'movingObjects_2D.obj movingGate_2D.obj \'
      write(22,FMT1)TAB,'movingPaddle_2D.obj movingWedge_2D.obj \'
      write(22,FMT1)TAB,'updateNormals_2D.obj vorticity_2D.obj\'
      write(22,FMT1)TAB,'periodicityCorrection_2D.obj \'


      !- Kernel Corrections
      if (i_kernelcorrection.eq.0) then
        write(22,FMT1)TAB,'ac_NONE_2D.obj \'
        write(22,FMT1)TAB,'kernel_correction_NC_2D.obj \'
      elseif (i_kernelcorrection.eq.2) then
        write(22,FMT1)TAB,'ac_KGC_2D.obj \'
        write(22,FMT1)TAB,'kernel_correction_KGC_2D.obj \'
        write(22,FMT1)TAB,'pre_self_KGC_2D.obj \'
        write(22,FMT1)TAB,'pre_celij_KGC_2D.obj \'
      elseif (i_kernelcorrection.eq.1) then
        write(22,FMT1)TAB,'ac_KC_2D.obj \'
        write(22,FMT1)TAB,'kernel_correction_KC_2D.obj \'
      endif

      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'ac_2D.obj \'
        write(22,FMT1)TAB,'poute_2D.obj \'
        write(22,FMT1)TAB,'gradients_calc_basic_2D.obj \'
      else
        write(22,FMT1)TAB,'ac_Conservative_2D.obj \'
        write(22,FMT1)TAB,'poute_Conservative_2D.obj \'
        write(22,FMT1)TAB,'gradients_calc_Conservative_2D.obj \'
      endif

      !- Boundary Conditions
      if(iBC.eq.1) then
        if(num_FB.eq.0)then
          write(22,FMT1)TAB,'monaghanBC_2D.obj \'
c          write(22,FMT1)TAB,'monaghanBC_AngMomCorrection_2D.obj \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,'self_BC_Monaghan_2D.obj \'
            write(22,FMT1)TAB,'celij_BC_Monaghan_2D.obj \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannConservative_2D.obj \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_RiemannNonConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_RiemannNonConservative_2D.obj \'
            endif
          endif
        else
          !- FLOATING BODIES
          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_2D.obj \'
c          write(22,FMT1)TAB,'monaghanBC_FloatingBodies_'
c     &                             //'AngMomCorrection_2D.obj \'
          if(iRiemannSolver.lt.1)then 
            write(22,FMT1)TAB,
     &                'self_BC_Monaghan_FloatingBodies_2D.obj \'
            write(22,FMT1)TAB,
     &               'celij_BC_Monaghan_FloatingBodies_2D.obj \'
          else
            if(iRiemannSolver.eq.1)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannConservative_2D.obj \'
            elseif(iRiemannSolver.eq.2)then 
              write(22,FMT1)TAB,
     &         'self_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.obj \'
              write(22,FMT1)TAB,
     &        'celij_BC_Monaghan_FloatingBodies_'
     &                               //'RiemannNonConservative_2D.obj \'
            endif
          endif
        endif
      else !(iBC.eq.2) 
        if(num_FB.eq.0)then
         write(22,FMT1)TAB,
     &        'self_BC_Dalrymple_2D.obj \'
         write(22,FMT1)TAB,
     &        'celij_BC_Dalrymple_2D.obj \'
        else
         write(22,FMT1)TAB,
     &        'self_BC_Dalrymple_FloatingBodies_2D.obj \'
         write(22,FMT1)TAB,
     &        'celij_BC_Dalrymple_FloatingBodies_2D.obj \'
        endif
      endif


      !- Rigid Body Motion 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'rigid_body_motion_2D.obj \'
      else
        write(22,FMT1)TAB,'rigid_body_motion_Conservative_2D.obj \'
      endif

      !- Variable Time Step 
      if(i_ConservativeFormulation.eq.0)then
        write(22,FMT1)TAB,'variable_time_step_2D.obj \'
      else
        write(22,FMT1)TAB,'variable_time_step_Conservative_2D.obj \'
      endif

      !- Approximate Riemann Solver
      if(iRiemannSolver.gt.0)then 
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'NonConservative_2D.obj \'
        else
          write(22,FMT1)TAB,'approx_RiemannSolver_'
     &                           //'Conservative_2D.obj \'
        endif
        i_limiterChoice = 1  !Temporary Default Value
        if(i_limiterChoice.eq.1)then
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.obj \'
        else
          write(22,FMT1)TAB,'limiter_betaMinMod_2D.obj \'
        endif
      endif

      !- Viscosity
      if (i_viscos.eq.1) then
          write(22,FMT1)TAB,'viscosity_artificial_2D.obj \'
      else if (i_viscos.eq.2) then
          write(22,FMT1)TAB,'viscosity_laminar_2D.obj \'
      else
          write(22,FMT1)TAB,'viscosity_laminar+SPS_2D.obj \'
      endif
		
      if(i_viscos.eq.3) then
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_SPS_2D.obj \'
        else
          print*
          print*,'Option not yet available: '
          print*,'SPS model in Conservative formulation'
          stop
          write(22,FMT1)TAB,'correct_SPS_Conservative_2D.obj \'
        endif
      else
        if(i_ConservativeFormulation.eq.0)then
          write(22,FMT1)TAB,'correct_2D.obj \'
        else
          write(22,FMT1)TAB,'correct_Conservative_2D.obj \'
        endif
      endif
        

c	    Kernels
c	    (1) Gaussian; (2)Quadratic;  (3)Cubic  ; (5) Wendland Quintic
      if (i_kernel.eq.1) then
             write(22,FMT1)TAB,'kernel_gaussian_2D.obj \'
      else if (i_kernel.eq.2) then 
             write(22,FMT1)TAB,'kernel_quadratic_2D.obj \'
      else if (i_kernel.eq.3) then 
             write(22,FMT1)TAB,'kernel_cubic_2D.obj \'
      else if (i_kernel.eq.5) then
	     write(22,FMT1)TAB,'kernel_wendland5_2D.obj \'
      endif

c         Equation of State
      if (i_EoS.eq.1) then
             write(22,FMT1)TAB,'EoS_Tait_2D.obj \'
      else if (i_EoS.eq.2) then
             write(22,FMT1)TAB,'EoS_IdealGas_2D.obj \'
      else if (i_EoS.eq.3) then
             write(22,FMT1)TAB,'EoS_Morris_2D.obj \'
      else
             write(22,FMT1)TAB,'EoS_Poisson_2D.obj \'
      endif

c         Density Reinitialisation
      if (i_densityFilter.eq.1) then
         write(22,FMT1)TAB,'densityFilter_Shepard_2D.obj \'
         write(22,FMT1)TAB,'ac_Shepard_2D.obj \'
         write(22,FMT1)TAB,'pre_celij_Shepard_2D.obj \'
         write(22,FMT1)TAB,'pre_self_Shepard_2D.obj \'
      else if (i_densityFilter.eq.2) then
         write(22,FMT1)TAB,'densityFilter_MLS_2D.obj \'
         write(22,FMT1)TAB,'ac_MLS_2D.obj \'
         write(22,FMT1)TAB,'LU_decomposition_2D.obj \'
         write(22,FMT1)TAB,'pre_celij_MLS_2D.obj \'
         write(22,FMT1)TAB,'pre_self_MLS_2D.obj \'
      else   
         write(22,FMT1)TAB,'densityFilter_NONE_2D.obj \'
      endif

c         Algorithm 
      if(i_ConservativeFormulation.eq.0)then
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_2D.obj '
        elseif (i_algorithm.eq.2) then   !Verlet

          write(22,FMT1)TAB,'step_verlet_2D.obj '
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_2D.obj '
        else

          write(22,FMT1)TAB,'step_beeman_2D.obj '
        endif
      elseif(i_ConservativeFormulation.eq.1)then   !ConservativeFormulation
        if (i_algorithm.eq.1) then       !Predictor Corrector
          write(22,FMT1)TAB,'step_predictor_corrector_'
     &                    //'Conservative_2D.obj '
        elseif (i_algorithm.eq.2) then   !Verlet
          write(22,FMT1)TAB,'step_verlet_2D.obj '
          print*
          print*,'Not a valid Option : '
          print*,'step_verlet in Conservative formulation'
          stop
        elseif (i_algorithm.eq.3) then   !Symplectic
          write(22,FMT1)TAB,'step_symplectic_Conservative_2D.obj '
        else
          write(22,FMT1)TAB,'step_beeman_2D.obj '
          print*
          print*,'Not a valid Option : '
          print*,'step_beeman in Conservative formulation'
          stop
        endif
      else
        print*
        print*,'Not a valid Option : '
        print*,'SPS model in Conservative formulation'
        stop
      endif
      write(22,FMT)'.f.obj:'
      write(22,FMT1)TAB,'ifort $(OPTIONS) $(COPTIONS) /O3 /c $<'
      write(22,FMT)      
      write(22,FMT)'SPHYSICS_2D.exe: $(OBJFILES)'
      write(22,FMT1)TAB,'xilink /OUT:$@ $(OPTIONS) $(OBJFILES)'

      write(22,FMT)      
      write(22,FMT)'clean:'
      write(22,FMT1)TAB,'del *.mod *.obj' !arno
       
      close(22)


      return
      end
	         
	   
	   
       
C-------SUBROUTINE TO ENSURE THAT PARTICLES ARE NOT TOO
C       CLOSE TO EACH OTHER
c
c     For Repulsive Boundary Force BC (iBC.eq.1), Boundary Particles
c     can be on top of each other
c     However, for Dalrymple BC (iBC.eq.2) particles ontop of each
c     will cause the code to crash
c
c
      subroutine position_check(dx,dz)
      include 'common.gen2D'
     
      dx_thresh=0.10*dx
      dz_thresh=0.10*dz
      
      
      !- Generate checking grid -
      one_over_2h = 1.0/(2.0*h)
      
      xGC_min = 0.0 - 4.0*h
      xGC_max = vlx
      yGC_min = 0.0 - 4.0*h
      yGC_max = vly
      zGC_min = 0.0 - 4.0*h
      zGC_max = vlz      
       
      ncx = int( (xGC_max-xGC_min)*one_over_2h ) + 1
      ncz = int( (zGC_max-zGC_min)*one_over_2h ) + 1
      
      nsheet  = ncx*ncz ! Number of cells in a XY sheet
      nct     = nsheet
      nc = 0
                   
      if(nct .ge. nct_max) then
        print*
        print*,'ERROR in position_check'
        print*,'No. of cells exceeds nct_max in common.gen2D'
        print*,'nct ',nct
        print*,'nct_max   ',nct_max
        print*
        stop
      endif 
                   
      !-- New 2-D Call numbers --
      ncall1  = -ncx+1               !South-East
      ncall2  =  ncx+1               !North-East
      ncall3  = -nsheet+1            !East       & Down
      ncall4  =  nsheet+1            !East       & Up
      ncall5  =  nsheet+ncx          !North      & Up
      ncall6  = -nsheet+ncx +1       !North-East & Down
      ncall7  =  ncall5+1            !North-East & Up
      ncall8  =  nsheet-ncx          !South      & Up
      ncall9  = -nsheet-ncx+1        !South-East & Down
      ncall10 =  ncall8+1            !South-East & Up
      
      print*
      print*,'-- Checking particles initial Positions --'
      print*,'dx, dz ',dx, dz
      print*,'dx_thresh ',dx_thresh
      print*,'dz_thresh ',dz_thresh
      print*,'ncx, ncz ',ncx,ncz

      call divide_dr_check(1,np)
      call ac_dr_check(dx_thresh,dz_thresh)

      print*,'--> Particle Positions Okay'
      print*

      return
      end
c     ____
  

c     ____________________ SUBROUTINE divide_dr_check

      
      subroutine divide_dr_check(n_start,n_end)
      include 'common.gen2D'
      
       
       do k=n_start,n_end
              
           dxp = xp(k) - xGC_min
           dzp = zp(k) - zGC_min
       
           icell = int( dxp * one_over_2h ) + 1
           kcell = int( dzp * one_over_2h ) + 1
       
           ! ii is the linear cell position in the matrix of 2h cells
           ii    = icell + (kcell - 1)*ncx 

          if(ii.lt.1)then
           print*
           print*,'ERROR in divide_2D.f'
           print*,'ii.lt.1'
           print*,'n_start,n_end ',n_start,n_end
           print*,'k',k
           print*,'xp(k) , zp(k)  ',xp(k), zp(k)
           print*,'xb_min, zb_min ',xb_min, zb_min
           print*,'dxp, dzp ',dxp, dzp
           print*,'one_over_2h  ',one_over_2h
           print*,'icell, kcell ',icell, kcell
           print*,'Box ii ',ii
           print*,'nc(ii)',nc(ii)
           stop
          end if
                                                        
           ! nc is the number of particles in cell ii
           nc(ii) = nc(ii)+1       
                                      
          if(nc(ii).gt.nplink_max)then
           print*
           print*,'ERROR in divide_dr_check'
           print*,'nc(ii) >= nplink_max'
           print*,'nc(ii)',nc(ii)
           print*,'nplink_max',nplink_max
           print*,'n_start,n_end ',n_start,n_end
           print*,'k',k
           print*,'xp(k), zp(k)',xp(k), zp(k)
           print*,'icell, kcell ',icell, kcell
           print*,'Box ii ',ii
           print*,'xb_min, zb_min ',xb_min, zb_min
           print*,'dxp, dzp ',dxp, dzp
           print*,'one_over_2h  ',one_over_2h
           stop
          end if
          
          ibox(ii,nc(ii))=k  !Tells us that particle with array location k (i.e. xp(k) )  
                                         !is in box ii which, so far, contains "nc(ii,mkind)" particles
           
       enddo 
       
       return
       end
c     ____
  

c     ____________________ SUBROUTINE ac_dr_check

      subroutine ac_dr_check(dx_thresh,dz_thresh)
      include 'common.gen2D'
      

      ncx_start = 1
      ncx_finish = ncx
c     ncy_start = 1
c     ncy_finish = ncy
      ncz_start = 1
      ncz_finish = ncz
      
       do lz=ncz_start,ncz_finish
       do lx=ncx_start,ncx_finish
       
       
         j1 = lx + (lz-1)*ncx

         if(nc(j1).gt.0) then
c                            ! if the cell is not empty, then
c                            ! loop over it and over neighboring
c                            ! cells

c         !-- Cells in the same XY sheet --
          lx2 = lx+1
          lz2 = lz+1
          
          !- EAST -
          if(lx2.le.ncx)then
            call celij(j1,j1+1,dx_thresh,dz_thresh)     !East
          end if

          !- NORTH -
          if(lz2.le.ncz)then
            call celij(j1,j1+ncx,dx_thresh,dz_thresh)     !North
          !- NORTH-EAST -  
            if(lx2.le.ncx)then 
            call celij(j1,j1+ncall2,dx_thresh,dz_thresh)     !North-East
            endif
          endif

          !- SOUTH-EAST -
          if(lx2.le.ncx.and.lz2.gt.2)then
            call celij(j1,j1+ncall1,dx_thresh,dz_thresh)     !South-East
          endif

          !-- PERIODIC CALLS ---
               
          !-- Cells in the same XZ sheet --
          !- X-Dir Periodic -
          if(i_periodicOBs(1).eq.1.and.lx.eq.ncx)then
            call celij(j1,j1+1-ncx,dx_thresh,dz_thresh)                 !East        Periodic X-Dir only
            if(lz2.le.ncz)then
              call celij(j1,j1+ncall2-ncx,dx_thresh,dz_thresh)          !North-East  Periodic X-Dir only
            endif
            if(lz2.gt.2)then
              call celij(j1,j1+ncall1-ncx,dx_thresh,dz_thresh)          !South-East  Periodic X-Dir only
            endif
            if(i_periodicOBs(3).eq.1)then
              if(lz.eq.ncz)then
                call celij(j1,j1+ncall2-ncx-nsheet,dx_thresh,dz_thresh) !North-East  Periodic X-Dir & Z-Dir
              else if(lz.eq.1)then
                call celij(j1,j1+ncall1-ncx+nsheet,dx_thresh,dz_thresh) !South-East  Periodic X-Dir & Z-Dir
              endif
            endif
          endif
          !- Z-Dir Periodic Only -
          if(i_periodicOBs(3).eq.1)then 
            if(lz.eq.ncz)then
              call celij(j1,j1+ncx-nsheet,dx_thresh,dz_thresh)          !North       Periodic Z-Dir only
              if(lx2.le.ncx)then
                call celij(j1,j1+ncall2-nsheet,dx_thresh,dz_thresh)     !North-East  Periodic Z-Dir only
              end if
            endif
            if(lz2.eq.2.and.lx2.le.ncx)then
              call celij(j1,j1+ncall1+nsheet,dx_thresh,dz_thresh)       !South-East  Periodic Z-Dir only
            end if
          end if

        endif  !End of: if(nc(j1,kind_p1).gt.0)then

      enddo
      enddo

       
       !do j1=1,nct
       do j1x = ncx_start, ncx_finish
       do j1z = ncz_start, ncz_finish
       
         j1 = j1x + (j1z-1)*ncx
           if(nc(j1).gt.0)
     +        call self(j1,dx_thresh,dz_thresh)
       enddo
       enddo
      
      end subroutine
c     ____
  

c     ____________________ SUBROUTINE celij_dr_check

      subroutine celij(j1,j2,dx_thresh,dz_thresh)
      include 'common.gen2D'
      
      
      do ii=1,nc(j1)
        i = ibox(j1,ii)
         
        do jj=1,nc(j2)
          j = ibox(j2,jj)

          drx = xp(i) - xp(j)
c         dry =
          drz = zp(i) - zp(j)
      
          delta_x=abs(drx)
c         delta_y=          
          delta_z=abs(drz)
          if(delta_x.lt.dx_thresh.and.
     &       delta_z.lt.dz_thresh) then
          
            print*
            write(*,*) 'WARNING: The following Particles are closer'
            write(*,*) 'than dx_thresh or dz_thresh', i,j
            print*,'xp(i), zp(i) ',xp(i), zp(i)
            print*,'xp(j), zp(j) ',xp(i), zp(j)
            if(delta_x.lt.0.1*dx_thresh.and.
     &         delta_z.lt.0.1*dz_thresh)then
     
              if((iBC.eq.1.and.(i.gt.nb.or.j.gt.nb)).or.iBC.eq.2)then
                print*,'Particles are too close and will cause code'
                print*,'code to crash'
c                stop
              end if   
                         
            end if
          endif
      
        enddo
      enddo
      
      end subroutine
      
c     ____
  
c     ____________________ SUBROUTINE self_dr_check

      subroutine self(j1,dx_thresh,dz_thresh)
      include 'common.gen2D'
      
      
      jj_start = 1
      do ii=1,nc(j1)
        i = ibox(j1,ii)
         
        jj_start = jj_start + 1
         
        do jj=jj_start,nc(j1)
          j = ibox(j1,jj)

          drx = xp(i) - xp(j)
c         dry =
          drz = zp(i) - zp(j)
      
          delta_x=abs(drx)
c         delta_y=          
          delta_z=abs(drz)
          if(delta_x.lt.dx_thresh.and.
     &       delta_z.lt.dz_thresh) then
          
            print*
            write(*,*) 'WARNING: The following Particles are closer'
            write(*,*) 'than dx_thresh or dz_thresh', i,j
            print*,'xp(i), zp(i) ',xp(i), zp(i)
            print*,'xp(j), zp(j) ',xp(i), zp(j)
            if(delta_x.lt.0.1*dx_thresh.and.
     &         delta_z.lt.0.1*dz_thresh)then
     
              if((iBC.eq.1.and.(i.gt.nb.or.j.gt.nb)).or.iBC.eq.2)then
                print*,'Particles are too close and will cause code'
                print*,'code to crash'
c                stop
              end if              
            end if
          endif
      
        enddo
      enddo
      
      end subroutine
      

C-------SUBROUTINE TO Choose precision of SPHYSICS variables:
c       xp, yp, zp, drx, dry, drz
c
c

      subroutine precisionWrite(i_compile_opt)
      include 'common.gen2D'
      
      print*
      print*,'Enter Desired Precision of XYZ variables'
      print*,'1 = Single Precision '
      print*,'2 = Double Precision '
      read(*,*)  i_Precision
      write(*,*) i_Precision
         
      if(i_Precision.eq.1)then
        kind_XZ = kind(0.0)
        print*,'Single Precision kind number is ',kind(0.0)
      elseif(i_Precision.eq.2)then
        kind_XZ = kind(0D0)
        print*,'Double Precision kind number is ',kind(0D0)
      else
        print*,'Invalid Choice for XYZ Precision'
        stop
      endif  
       
      print*,'Opening source\SPHYSICS2D\precision_kind.2D, '
      if(i_compile_opt.ge.3)then   !Windows Choice
        open(unit=19,file='..\..\source\SPHYSICS2D\precision_kind.2D')
      else
        open(unit=19,file='../../source/SPHYSICS2D/precision_kind.2D')
      endif
      
      write(19,190)'      parameter(kind_XZ=',kind_XZ,')'
190   format(a24,i2,a1)

      close(19)
      print*
      
      return
      end