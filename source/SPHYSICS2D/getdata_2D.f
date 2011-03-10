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

      subroutine getdata
      character*10  ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10,ch11,ch12

      include 'common.2D'
      

c

      character supp_ini*4,name_ini*40,normals_ini*40
      character (LEN=40) :: paddle_fileName, MovingObjectsFile 
c
      open(11,file='INDAT')
 
      read(11,*)i_kernel
      read(11,*)i_algorithm
      read(11,*)i_densityFilter  ! Using a density filter
      read(11,*)i_viscos
      read(11,*)iBC              !BC
      read(11,*)i_periodicOBs(1)
      read(11,*)i_periodicOBs(2)
      read(11,*)i_periodicOBs(3)
      read(11,*)lattice
      read(11,*)i_EoS
      read(11,*)h_SWL 
      read(11,*)B
      read(11,*)gamma
      read(11,*)coef
      read(11,*)eps   ! epsilon XSPH
      read(11,*)rho0
      read(11,*)viscos_val
      read(11,*)visc_wall
      read(11,*)vlx    !Dimensions of container box
      read(11,*)vly    ! vly=0 for 2D
      read(11,*)vlz
      read(11,*)dx0
      read(11,*)dy0     ! dy=0  for 2D
      read(11,*)dz0
      read(11,*)h
      read(11,*)np	
      read(11,*)nb
      read(11,*)nbf
      read(11,*)ivar_dt
      read(11,*)dt
      read(11,*)tmax     ! End of the run(seconds)
      read(11,*)out      ! output every out seconds
      read(11,*)trec_ini ! Initial output 
      read(11,*)dtrec_det !Detailed recording
      read(11,*)t_sta_det
      read(11,*)t_end_det
      read(11,*)i_restartRun
      read(11,*)CFL_number
      read(11,*)TE0
      read(11,*)i_kernelcorrection
      read(11,*)iRiemannSolver 
      read(11,*)iTVD
      read(11,*)beta_lim
      read(11,*)i_vort
      read(11,*)ndt_VerletPerform
      read(11,*)ndt_FilterPerform
      read(11,*)ndt_DBCPerform
 
          
       write(80,*) ' '
       write(80,*)'Periodicity '
       write(80,*)'i_periodicOBs'
       write(80,*)'X-Direction: ',i_periodicOBs(1)
       write(80,*)'Y-Direction: ',i_periodicOBs(2)
       write(80,*)'Z-Direction: ',i_periodicOBs(3)
       write(80,*) ' '
       !- Screen printout -
       print*
       print*,'Periodicity '
       print*,'i_periodicOBs'
       print*,'X-Direction: ',i_periodicOBs(1)
       print*,'Y-Direction: ',i_periodicOBs(2)
       print*,'Z-Direction: ',i_periodicOBs(3)
       print*
       


       if (np .ge. npar ) then
         write(80,*) 'Number of particles exceeds npar in common!'
         write(*,*) ' '
         write(*,*) 'Number of particles exceeds npar in common!'
         stop
       endif
       if (nb .ge. nb_max ) then
         write(80,*) 'Number of Boundary particles exceeds'
         write(80,*) 'nb_max in common!'
         write(*,*) 'Number of Boundary particles exceeds'
         write(*,*) 'nb_max in common!'
         stop
       endif
       write(80,*) 'np, nb, nbf, ivar_dt, dt ',np,nb,nbf,ivar_dt,dt
       write(80,*) 'h = ',h
       write(80,*) '(dx0,dy0,dz0)=',dx0,dy0,dz0
       write(80,*) B,viscos_val,gamma, rho0
       write(80,*) 'Still Water Depth h_SWL = ',h_SWL
       write(80,*) 'CFL_number = ',CFL_number
       if(iTVD.eq.1)write(80,*) 'Using TVD, beta_lim =  ',beta_lim
       !-- screen printout --
       write(*,*) 'np, nb, nbf, ivar_dt, dt ',np,nb,nbf,ivar_dt,dt
       write(*,*) 'h = ',h
       write(*,*) '(dx0,dy0,dz0)=',dx0,dy0,dz0
       write(*,*) B,viscos_val,gamma, rho0
       write(*,*) 'Still Water Depth h_SWL = ',h_SWL
       write(*,*) 'CFL_number = ',CFL_number
       if(iTVD.eq.1)write(*,*) 'Using TVD, beta_lim =  ',beta_lim
       
       
       close(11)


       nbfp1=nbf+1
       nbp1=nb+1
       
c     speed of sound at the reference density 

       if(i_EoS.eq.1)then
         cs0 = sqrt(gamma*B/rho0)
         write(80,*)'cs0=sqrt(gamma*B/rho0)= ',cs0
         write(*,*) 'cs0=sqrt(gamma*B/rho0)= ',cs0
         i_gamma = nint(gamma)
         write(80,*)'Integer gamma, i_gamma = nint(gamma) = ',i_gamma
         write(*,*) 'Integer gamma, i_gamma = nint(gamma) = ',i_gamma
       endif

       ddt_c=dt   ! To be used with variable time step
       ddt_p=dt   ! and predictor corrector
       
       
       if(i_algorithm.eq.2)then
         nstep_verlet = 0  !Initialisation of verlet counter
       end if

       if(iBC.eq.2)then	!Hughes & Grahams (2009) correction 
         nstep_DBC=0       !Initialisation of DBC correction counter
         write(80,*)'ndt_DBCPerform ',ndt_DBCPerform
         write(*,*)'ndt_DBCPerform ',ndt_DBCPerform
       else
         nstep_DBC=-1      !when iBC=1
       endif
       
c      -- Parameters for repulsive force BC --       
       if(lattice.eq.1)then
         vnorm_mass = 1.0
       else if(lattice.eq.2)then
         vnorm_mass = 0.5
       end if
       write(80,*)'lattice    ',lattice
       write(80,*)'vnorm_mass ',vnorm_mass
       write(*,*) 'lattice    ',lattice
       write(*,*) 'vnorm_mass ',vnorm_mass
       if(iBC.eq.1)then
         beta_coef2 = (coef/10.)**2
         beta_coef  = (coef/10.)
         one_over_beta_coef2 = 1.0/beta_coef2
         one_over_beta_coef  = 1.0/beta_coef
       else
         one_over_beta_coef2 = 0.
         one_over_beta_coef  = 0.
       end if

       if(i_densityFilter.gt.0)then
        write(80,*)'Initialising ndt_FilterPerform ',ndt_FilterPerform
        write(*,*) 'Initialising ndt_FilterPerform ',ndt_FilterPerform
       endif

C      Initializing thermal energy
       do i=1,np
         TEp(i)=0.
       enddo    
      
C      Initial image to RUN
     
c      i_restartRun > 1 is used for Checkpointing = repetitive restarting of code (for clusters)
c        i_restartRun = 0   :   Start new run, once only
c        i_restartRun = 1   : reStart old run, once only
c        i_restartRun = 2   :   Start new run, with repetitive restarts (Checkpointing)
c        i_restartRun = 3   : reStart old run, with repetitive restarts (Checkpointing)
 
       if(i_restartRun.ge.1) then   !if(i_restartRun.eq.1) then  

         i_RESTART_exists = 1
         open(44,file='RESTART',err=45)
         read(44,100,err=45)itime,time,ngrab,dt
         close(44)

         i_RESTART_exists = i_RESTART_exists + 1
45       i_RESTART_exists = i_RESTART_exists - 1


         if(i_RESTART_exists.eq.1.and.itime.gt.0)then
           write(supp_ini,'(i4.4)') ngrab
           name_ini='PART_'//supp_ini
           write(*,*) name_ini
           normals_ini='NORMALS.init'   !Used for compatibility with parallel Code
           MovingObjectsFile='Floating_Bodies.RESTART'
         else
           itime = 0
           ngrab = 0
           grab_P = 0.0
           grab_E = 0.0
           name_ini='IPART'
           write(*,*) name_ini
           normals_ini='NORMALS.init'
           MovingObjectsFile='Floating_Bodies.txt'
         endif

       else
         itime = 0
         ngrab = 0
         grab_P = 0.0
         grab_E = 0.0
         name_ini='IPART'
         write(*,*) name_ini
         normals_ini='NORMALS.init'
         MovingObjectsFile='Floating_Bodies.txt'

       endif
100    format(i8,e16.8,i8,e16.8)
       write(80,*)'i_restartRun ',i_restartRun 
       write(80,*)'Reading in ', name_ini,normals_ini
       write(80,*)'Using dt ', dt  !BDR
       write(*,*) ' '
       write(*,*)'i_restartRun ',i_restartRun 
       write(*,*)'Reading in ', name_ini,normals_ini
       write(*,*)'Using dt ', dt  !BDR


       
c     --- Load Particle Position Data ---  
       
       open(13,file=name_ini)

196    format(7e16.8)

       do i=1,np
         read(13,196) xp(i),zp(i),up(i),wp(i),rhop(i),p(i),pm(i)
                                                                
         call equation_of_state(rhop(i),TEp(i),p(i),cs(i))

         if(i_algorithm.eq.2) then ! To be used with Verlet algorithm
           xo(i)=xp(i)            
           zo(i)=zp(i)
           rhoo(i)=rhop(i)
           uo(i)=up(i)
           wo(i)=wp(i)
           pr(i) = p(i)/(rhop(i)*rhop(i))
           TEo(i)=TEp(i)
         elseif (i_algorithm.eq.4) then  !to be used with Beeman algorithm
           udot(i)=0.0   !Check it for hot starts !!!!!!!            
           wdot(i)=0.0
           udoto(i)=0.0
           wdoto(i)=0.0
           udotm1(i)=0.0
           wdotm1(i)=0.0
           rdot(i)=0.0
           rdoto(i)=0.0
           rdotm1(i)=0.0

         endif
       
         iflag(i)=1 ! To detect particles inside the container

         pVol(i)      = pm(i)/rhop(i)
         
         if(iRiemannSolver.eq.1)then
           Volrho(i)    = pm(i)
           Volrhou(i)   = pVol(i)*rhop(i)*up(i)
           Volrhow(i)   = pVol(i)*rhop(i)*wp(i)
           TotalE(i)    = rhop(i)*(TEp(i)
     &                  + 0.5*(up(i)*up(i)+wp(i)*wp(i)))
           VolTotalE(i)  = pVol(i)*TotalE(i)
         endif    
            
       enddo
       close(13) 
c      ----------------------------------
       
c     ---  Open NORMALS Files ---
      if(iBC.eq.1) then
        open(15,file=normals_ini)       
        do i=1,nb
          read(15,*) 
     &    xnb(i),znb(i),
     &    iBP_Pointer_Info(i,1),       !Absolute index BP     
     &    iBP_Pointer_Info(i,2),       !Rank of BP (default=0, reserved for MPI)       
     &    iBP_Pointer_Info(i,3),       !Absolute index of i-1 neighbour BP      
     &    iBP_Pointer_Info(i,4),       !Absolute index of i+1 neighbour BP      
     &    BP_xz_Data(i,1),
     &    BP_xz_Data(i,2)
        enddo
        nstart=nbp1
        close(15)
        call updateNormals_2D(1,nb)
        !- Angular Momentum Correction -
        two_alpha_h = 2.0*1.3
        rNum_h = 3.0*h
        one_over_rNum_h = 1.0/rNum_h 
        one_over_r_ij_dot_n = 1.3 /h
      else
        nstart=1
      endif
      nstart_minus1 = nstart - 1
c     -----------------------------
              
              
c     ------ Load Moving Object Data Files ---------
c
c     - (Note files will be empty for no objects) -
              

c     --- The gate ----------------
       open(77,file='gate')
       read(77,*) iopt_gate       
       if(iopt_gate.eq.1) then
         read(77,*) ngate_ini
         read(77,*) ngate_end
         read(77,*) VXgate,VZgate
         read(77,*) tgate
       else
         ngate_ini=1  ! If ngate_ini > ngate_end the loop 
         ngate_end=0  ! in AC doesn't work 
       endif       
       close(77)

       !-- The wavemaker --
       open(66,file='wavemaker')
       read(66,*) iopt_wavemaker       
       if(iopt_wavemaker.eq.1)then
         read(66,*) i_paddleType
         read(66,*) nwavemaker_ini
         read(66,*) nwavemaker_end
         read(66,*) X_PaddleCentre
         read(66,*) X_PaddleStart
         read(66,*) paddle_SWL
         read(66,*) flap_length
         read(66,*) stroke
         read(66,*) twavemaker
         read(66,*) paddle_fileName
         read(66,*) Nfreq
       
      
         do n=1,Nfreq
           read(66,*) A_wavemaker(n)
           read(66,*) Period(n)
           read(66,*) phase(n)
           read(66,*) twinitial(n)
         enddo
         
         do i=1,nb
           xp_ini(i)=xp(i)  ! X initial of two layers of wavemaker
         end do
       
         A_wavemaker(Nfreq+1)=A_wavemaker(Nfreq)
         Period(Nfreq+1)=Period(Nfreq)
         phase(Nfreq+1)=0.
         twinitial(Nfreq+1)=100*tmax
       
         ind_wm=0
         ind1_wm=1
         
         if(i_paddleType.eq.3)then  !Paddle has prescribed motion
           print*
           print*,'Opening prescribed paddle file = ',paddle_fileName
           open(67,file=paddle_fileName)
           n_paddleData = 0
           i_read = 0
           do while(i_read.eq.0)
             n_paddleData = n_paddleData + 1
             if(n_paddleData.gt.20000+1)then
               print*,'Too many timesteps in prescribed paddle motion'
               stop
             endif
             read(67,*,END = 68)time_paddle(n_paddleData),
     &                x_paddle(n_paddleData),u_paddle(n_paddleData)
             if(i_read.eq.0)then
               i_read = i_read - 1
             endif
68           i_read = i_read + 1
           enddo
           close(67)
           n_paddleData = n_paddleData - 1
           print*
           print*,'Sample print-out of prescribed paddle motion'
           print*,'  n_paddleData = ',n_paddleData
           print*,'  i,  time_paddle(i),    x_paddle(i),   u_paddle(i)'
           do i=1,10
             print*,i,time_paddle(i),x_paddle(i),u_paddle(i)
           enddo
           print*
         endif
              
       else
       
         nwavemaker_ini=1  ! If nwavemaker_ini > nwavemaker_end the loop 
         nwavemaker_end=0  ! in AC doesn't work 
       
       endif       
       close(66) 
c     -----------------------------

c      ---  Open Raichlen Wedge Input Files --
       open(77,file='Tsunami_Landslide.txt')
       read(77,*) iopt_RaichlenWedge
       if(iopt_RaichlenWedge.eq.1)then
         read(77,*)bslope
         write(80,*)' '
         write(80,*)'-- Opening Raichlen Wedge Input Files --'
         write(80,*)' '
         write(*,*) ' '
         write(*,*) '-- Opening Raichlen Wedge Input Files --'
         write(*,*) ' '
         open(unit=78,file='Benchmark_4_run30.txt')
         read(78,*) i_RBend
         write(80,*)'Number of time measurements, i_RBend ',i_RBend
         write(*,*) 'Number of time measurements, i_RBend ',i_RBend
         do i = 1,i_RBend
           read(78,*)RBtime(i),RBxpos(i)
           RBxpos(i) = RBxpos(i)*0.01   !conversion from cm to m
         end do
         close(78)
         do i = nbf+1,nb
           xp_Wedge_initial(i-nbf) = xp(i)
           zp_Wedge_initial(i-nbf) = zp(i)

         end do
         write(80,*)'Wedge info loaded'
         write(80,*)' '
         write(*,*) 'Wedge info loaded'
         write(*,*) ' '
       end if
       close(77) 
c      ---------------------------------------
       
c      ---  Open Floating Body Input Files --
       open(77,file=MovingObjectsFile)  !open(77,file='Floating_Bodies.txt')
       read(77,*) iopt_FloatingBodies
       i_FB_Pointer_Info = 0
       if(iopt_FloatingBodies.eq.1)then
         read(77,*)nbfm
         print*
         print*,'-- Opening Floating Body Input Files --'
         i_FB_finish = 0
         i_ini = nbfm + 1
         num_FB_counter = 0
         do while(i_FB_finish.eq.0)
           read(77,*,END = 76)i_num_FB
           if(i_num_FB.gt.num_FB_max)then
             print*,'Number of floating Bodies exceeds max value'
             print*,'num_FB.gt.num_FB_max'
             print*,'Adjust num_FB_max in common.2D, num_FB_max = ',
     &          num_FB_max
             stop
           endif
           read(77,*)bigMass(i_num_FB)
           read(77,*)bigInertiaYY(i_num_FB)
c
c
           read(77,*)XcylinderDimension(i_num_FB),
     &               ZcylinderDimension(i_num_FB)
           read(77,*)cylinderDensity(i_num_FB)
           read(77,*)FB_SpecificWeight(i_num_FB)
           read(77,*)friction_coeff(i_num_FB)
           read(77,*)Box_XC(i_num_FB),Box_ZC(i_num_FB)
           read(77,*)bigU(i_num_FB),bigW(i_num_FB),
     &                      bigOmegaY(i_num_FB)
           read(77,*)nb_FB(i_num_FB)
           print*,'i_num_FB     ',i_num_FB
           print*,'bigMass      ',bigMass(i_num_FB)
           print*,'bigInertiaYY ',bigInertiaYY(i_num_FB)
           print*,'XcylinderDimension ',XcylinderDimension(i_num_FB)
           print*,'ZcylinderDimension ',ZcylinderDimension(i_num_FB)
           print*,'cylinderDensity    ',cylinderDensity(i_num_FB)
           print*,'SpecificWeight     ',FB_SpecificWeight(i_num_FB)
           print*,'friction_coeff     ',friction_coeff(i_num_FB)
           print*,'Initial Conditions:'
           print*,'Box_XC, Box_ZC    ',
     &             Box_XC(i_num_FB),Box_ZC(i_num_FB)
           print*,'bigU, bigW  ',
     &             bigU(i_num_FB), bigW(i_num_FB)
           print*,'OmegaY  ', 
     &             bigOmegaY(i_num_FB)
           print*,'nb_FB        ',nb_FB(i_num_FB)
           do i = i_ini,nb_FB(i_num_FB)
             i_FB_Pointer_Info(i) = i_num_FB
             xp_minus_R0(i) = xp(i) - Box_XC(i_num_FB)
             zp_minus_R0(i) = zp(i) - Box_ZC(i_num_FB)

           end do
c          - zero body forces -
           BodyForce_x      = 0.0
           BodyForce_z      = 0.0              

           X_nonFriction = 0.0
           Z_nonFriction = 0.0

           i_ini = nb_FB(i_num_FB) + 1
           num_FB_counter = num_FB_counter + 1
           u_parallel_max = cs0*1.0e-4  !Threshold sliding velocity to activate dynamic friction force
           one_over_u_parallel_max = 1.0 / u_parallel_max
           
           !Determine whether to exit loop
           if(i_FB_finish.eq.0)then
             i_FB_finish = i_FB_finish - 1
           endif
76           i_FB_finish = i_FB_finish + 1
           print*
           
         end do
         num_FB = num_FB_counter
         print*,'Number of Floating Bodies, num_FB ',num_FB
         print*,'nbfm, nstart ',nbfm,nstart
         print*,'--- Floating Body info loaded ---'
         print*
       end if
       close(77) 
       !stop
c      ---------------------------------------
       
       if(iopt_gate          .eq.1.or.
     &    iopt_wavemaker     .eq.1.or.
     &    iopt_RaichlenWedge .eq.1.or.
     &    iopt_FloatingBodies.eq.1) then
         iopt_movingObject = 1
       else
         iopt_movingObject = 0           
       endif
       write(80,*) ' '
       write(80,*) 'Moving Objects:iopt_movingObject ',iopt_movingObject
       write(80,*) ' '
       write(*,*) ' '
       write(*,*) 'Moving Objects:iopt_movingObject ',iopt_movingObject
       write(*,*) ' '

c     --- End of Loading Moving Object Files ---


       !--  Sub-Particle Scale (SPS) Turbulent Stress Terms --     
       !- Smag = (0.12*dp)**2  -
       dp=sqrt(dx0*dx0 + dz0*dz0)/2.0
       Smag = (0.12*dp)**2 
       twoThirds = 2.0/3.0
       Blin_const = twoThirds*0.0066*dp*dp       
       
C
       write(80,*)' All ',np,' initial points read successfully'
       write(*,*)' All ',np,' initial points read successfully'
       
       pi = 4.d0*datan(1.d0)
       twopi=2.*pi

C     The extent of all particles is the same but it could be different
c    (1= boundary-boundary; 3= fluid-fluid; 2 boundary-fluid)
 
       write(80,*) ' '
       write(80,*)' -- Kernel Information --'
       write(80,*)'Kernel Choice = ',i_kernel
       write(80,*) ' '
       
       !Different Kernel Smoothing Lengths

       two_h  = 2.0*h            
       four_h = 4.0*h
       six_h  = 6.0*h    
       h_over_100 = 0.01*h        
       h_over_1000 = 0.001*h        
 
c     2D Kernel
c.... coefficients for the spline Kernel W(r)
c.... coefficients for the derivative W_r (r,h) of the spline Kernel
c     (See Liu's Book 2003)	      
       
         h2=h*h
         fourh2=4.*h2

       
        if (i_kernel.eq.1) then !2D Gaussian 
c                    Wab=      a2*exp(-q2/gauss2)              gaussian kernel
c                              0
c                    W_r=      f1*q*exp(-q2/gauss2)   derivative of gaussian kernel
c                              0
           gauss=1.
           gauss2=gauss*gauss                           
           a1 = 1./pi
           a2 = a1/(h*h)
           aa = a1/(h*h*h)
           adh=a2
           f1 = -2.*aa/gauss2       
       
           write(80,*) 'Gaussian kernel'
           write(80,*) a2,f1
       
           deltap2=0.5
           Wdeltap=a2*exp(-deltap2/gauss2)
           od_Wdeltap=1./Wdeltap
           write(80,*) od_Wdeltap
       
       
        elseif (i_kernel.eq.2) then !2D Quadratic 
c               Wab=      b1*(0.25*q*q-q+1)       	quadratic spline-kernel
c               0
c               W_r=      e1*(0.5*q-1)   derivative of quadratic spline-kernel
c               0
           a1 = 10./(pi*7.)
           a2 = a1/h2
           aa = a1/(h2*h)
              
           b1 = 21.*a2/20.
           e1 = 21.*aa/20.
           adh=b1
           write(80,*) 'Quadratic Kernel',aa
           write(80,*) 'b1,e1 ',b1,e1

        elseif(i_kernel.eq.3)then !2D Cubic Spline
       
c             Wab=      a2*(1-1.5*q*q+0.75*q*q*q)
c                       a24*(2-q)*(2-q)*(2-q)           cubic spline-kernel
c                       0
c             W_r=      c1*q + d1*q*q
c                       c2*(2-q)**2   derivative of cubic spline-kernel
c                       0
       
           a1 = 10./(pi*7.)
           a2 = a1/h2
           aa = a1/(h2*h)
           adh=a2  
           a24 = 0.25*a2
       
           c1 = -3.*aa
           d1 = 9.*aa/4.
           c2 = -3.*aa/4.
       
           write(80,*) 'cubic spline',aa,a24 
           write(80,*) c1,a2,d1

c                 W(Delta_p): used to calculate fab to eliminate the
c                          Tensile Instability. We'll consider
c                 Delta_p=h/1.3 as suggested by Monaghan

           deltap=1./1.5    
           Wdeltap=a2*(1.-1.5*deltap*deltap+
     +        0.75*deltap*deltap*deltap)
           od_Wdeltap=1./Wdeltap
           write(80,*) od_Wdeltap
       
       elseif (i_kernel.eq.5) then  !2D Wendland Quintic
       
       
!		W=Awen*(1-qq/2)**4*(2*qq+1) 0<=r<=2*h
!		dw/dr=Bwen*qq*(1-qq/2)**3. 0<=r<=2*h
!		Awen(i)=7/(4*pi*h3)
!		Bwen(i)=-35/(4*pi*h4)

           Awen=0.557/(h2)
           Bwen=-2.7852/(h2*h)
           adh=Awen                            
           write(80,*) 'Wendland kernel'
           write(80,*) Awen,Bwen
       
           deltap=0.5
           Wdeltap=Awen*((1.0-deltap*0.5)**4.0)*(2.0*deltap+1.0)
           od_Wdeltap=1./Wdeltap
           write(80,*) od_Wdeltap
             
       else
           write(80,*) 'i_kernel=',i_kernel,
     +                 ' is not a valid option !!!'
           write(*,*) 'i_kernel=',i_kernel,
     +                ' is not a valid option !!!'
       endif

c
c ,,,, viscosity term eta2
c
       eta  = 0.1*h
       eta2 = eta*eta
       
c
c  ... determine xmin,xmax,ymin,ymax,zmin,zmax
c
       xmin = xp(1)
       xmax = xp(1)
       zmin = zp(1)
       zmax = zp(1)
      
       do i=2,nb
         if(xp(i).gt.xmax) xmax=xp(i)
         if(xp(i).lt.xmin) xmin=xp(i)
         if(zp(i).gt.zmax) zmax=zp(i)
         if(zp(i).lt.zmin) zmin=zp(i)
       enddo

       !- Limits of domain -   
       xmax_container = xmax 
       xmin_container = xmin 
       zmin_container = zmin
       
       !Periodic correction
       if(i_periodicOBs(1).eq.1)then
         xmax_container = xmax + vnorm_mass*dx0   !BDR
       endif
       xmax_container_double = dble(xmax_container)
       xmin_container_double = dble(xmin_container)

       write(80,*) ' '
       write(80,*)'xmin_container, xmax_container = ',
     &             xmin_container, xmax_container
       write(80,*)'zmin_container = ',zmin
       write(80,*) ' '
       write(*,*) ' '
       write(*,*)'xmin_container, xmax_container = ',
     &            xmin_container, xmax_container
       write(*,*)'zmin_container = ',zmin
       write(*,*) ' '
       
       do i=nbp1,np
         if(xp(i).gt.xmax) xmax=xp(i)
         if(xp(i).lt.xmin) xmin=xp(i)
         if(zp(i).gt.zmax) zmax=zp(i)
         if(zp(i).lt.zmin) zmin=zp(i)
       enddo
       
c      -- Add small border region to prevent error --
c      -- when calculating icell.. in divide       --
       if(i_periodicOBs(1).ne.1)then
         xmax = xmax + 0.1*h
         
       endif
       xmin = xmin - 0.1*h

       zmax = zmax + 0.1*h
       zmin = zmin - 0.1*h
       
       if(iopt_wavemaker.eq.1.and.i_paddleType.eq.2)then
         !xmin for piston-flap paddle extends further back
         !than initial xmin
         print*,'Correcting xmin for paddle '
         print*,'Current xmin   = ',xmin
         print*,'X_PaddleCentre = ',X_PaddleCentre
         print*,'stroke         = ',stroke
         print*
         xmin_check = X_PaddleCentre - 1.5*stroke
         if(xmin_check.lt.xmin)then
           xmin = xmin_check
         endif
       end if
       
       zmax_ini=zmax  ! used to keep moving obstacles 
                      ! inside the initial boundaries
       
       ! To place particles in a convenient place that leave the domain
       xtrash=2.*xmax
       ztrash=0.

c
c  ... determine number of cells
c
       one_over_h  = 1.0/h
       one_over_2h = 1.0/(2.0*h)
       
       ncx = int( (xmax-xmin)*one_over_2h ) + 1
       ncz = int( (zmax-zmin)*one_over_2h ) + 1

       
        ! 2b used in check_limits 
       zcontrol=zmin+ncz*2*h
       write(80,*)'zcontrol ',zcontrol,zmax_ini

       !- Periodic Correction -
       if(i_periodicOBs(1).eq.1)then
         ncx_min = int((xmax_container-xmin)*one_over_2h)
         if(ncx.gt.ncx_min)ncx = ncx_min
         write(80,*)'Corrected ncx for X-Periodicity'
         write(*,*) 'Corrected ncx for X-Periodicity'
         write(*,*) 'ncx, ncx_min', ncx, ncx_min
         write(*,*) 'xmin, xmax_container', xmin, xmax_container
         write(*,*) 'xmin_container_double, xmax_container_double ',
     &               xmin_container_double, xmax_container_double
       endif
c       if(i_periodicOBs(3).eq.1)then
c         ncz_min = int((zmax_container-zmin)*one_over_2h)
c         if(ncz.gt.ncz_min)ncz = ncz_min
c         write(80,*)'Corrected ncz for Z-Periodicity'
c         write(*,*) 'Corrected ncz for Z-Periodicity'
c         write(*,*) 'ncz, ncz_min'
c         write(*,*) 'zmin, zmax_container',zmin, zmax_container
c       endif

       nct = ncx*ncz
       
       if(nct.gt.nct_max)then
         write(80,*) ' '
         write(80,*)'ERROR in getdata.f'
         write(80,*)'No. of cells exceeds nct_max in common.2D'
         write(80,*)'nct ',nct
         write(80,*)'nct_max   ',nct_max
         write(80,*) ' '
         write(*,*) ' '
         write(*,*)'ERROR in getdata.f'
         write(*,*)'No. of cells exceeds nct_max in common.2D'
         write(*,*)'nct ',nct
         write(*,*) ' '
         stop
       endif
              
       ncz_ini=ncz
       nct_ini=nct
       
       write(80,*)'xmin, xmax = ',xmin,xmax
       write(80,*)'zmin, zmax = ',zmin,zmax
       
       write(80,*)'ncx/z    = ',ncx,ncz
       write(80,*)'nct= ',nct

c      
c     --- Periodic ncall numbers for X-direction periodicity ---
       ncall11 =     1 - ncx          !East X-Periodic
       ncall12 = ncx-1 + ncx    !North-West X-Periodic
       ncall13 = ncx+1 - ncx    !North-East X-Periodic

       return
       end
       
c     ___________________
