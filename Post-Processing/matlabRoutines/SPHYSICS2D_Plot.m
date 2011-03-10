% POST PROCESSING ROUTINE FOR SPHYSICS_2D. THIS PROGRAM GENERATES PARTICLE
% PLOTS OF THE RESULTS USING THE PART_ijkl,DT AND matlabin FILES AS INPUT.
% THIS PROGRAM SEARCHES FOR THE PART_ijkk, DT AND matlabin FILES IN THE 
% CURRENT MATLAB WORKING DIRECTORY. 

clear all
close all

%READ PARAMETERS FROM MATLABIN
load matlabin;
np=matlabin(1)
vlx=matlabin(2)
vly=matlabin(3)
vlz=matlabin(4)
out=matlabin(5)
nb=matlabin(6)
nbf=matlabin(7)

%IF THE SIMULATION IS COMPLETE LOAD THE TIME INFORMATION FROM FILE DT. IF
%SIMULATION IS ONGOING, JUST PLOT USER SPECIFIED FRAMES
irun=input('Is the run complete? (1=Yes,0=No)  ');
if(irun==1)
    %LOAD AND READ TIME,DT FROM FILE DT. THE FIRST HEADLINE IN THE FILE IS
    %SKIPPED
    [time,DT1,DT2,DT]=textread('DT','%f %f %f %f','headerlines',1);
    %THE NUMBER OF PARTFILES, Nframes, IS ONE LESS THAN THE LENGTH OF THE time VECTOR.
    %THIS IS BECAUSE THE LAST VALUE OF THE TIME VECTOR IS ASSOCIATED WITH EPART
    Nframes_ini=1;
    Nframes=length(time)-1;
else
    Nframes_ini=input('Start Number of PART files to be visualized= ');
    Nframes=input('Number of PART files to be visualized= ');
end

sframes=input('Save images in jpg format ? (1=Yes,0=No)  ');

%INITIALIZE FIGURE AND SET PROPERTIES
XF=figure;
set(XF,'DefaultAxesFontSize',16)
set(XF,'DefaultTextFontSize',16)

%SET NORMALIZED FIGURE POSITION AND SIZE TO MAKE IT CONSISTENT OVER ALL
%MONITORS
set(XF,'Units','normalized')
%set(gcf,'Position',[0.1,0.2,0.8,0.7])  %Dam Break - Case 1
set(gcf,'Position',[0.05,0.5,0.95,0.4])  %2-D Waves - Case 3 (modified)
set(gca,'Position',[0.07 0.14 0.90 0.75])
%set(gcf,'Position',[0.1,0.5,0.8,0.4])   %Original


%INITIALIZE PARTICLE SIZE, VARIABLE M FOR MOVIES, PLOT AXES AND
%POSITION OF TEXT
ParticleSize=12;
M=moviein(Nframes);
xmin=-0.25   %-0.25  %-0.01;
xmax=vlx   %2.25   %1.5*vlx;
zmin=-0.01  %-0.25  %-0.01;
zmax=1.5*vlz   %Original
zmax=1.25*vlz  %Dam Break
zmax=1.25*vlz   %Waves
i_geometry=input('Zoom? (1=Yes,0=No)  ');
if i_geometry == 1
  xmin = 16.3
  xmax = 17.9
  zmin = 0.3
  zmax =  0.7
end

% IF THE SIMULATION HAS FIXED BOUNDARIES CALL NO MOVING BOUNDARIES ELSE
% CALL MOVING BOUNDARIES. VARIABLE nbf IS THE NUMBER OF FIXED BOUNDARY
% PARTICLES
if(nbf<nb) 
    MovingObjects2D
else
    NoMovingObjects2D
end
