% POST PROCESSING ROUTINE FOR SPHYSICSgen_2D. THIS PROGRAM PLOTS THE INITIAL
% PARTICLE CONFIGURATION USING IPART. 
% THIS PROGRAM SEARCHES FOR THE IPART, DT AND matlabin FILES IN THE 
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

%INITIALIZE FIGURE AND SET PROPERTIES
XF=figure;
set(XF,'DefaultAxesFontSize',16)
set(XF,'DefaultTextFontSize',16)

%SET NORMALIZED FIGURE POSITION AND SIZE TO MAKE IT CONSISTENT OVER ALL
%MONITORS
set(XF,'Units','normalized')
%set(gcf,'Position',[0.1,0.2,0.8,0.7])  %Dam Break - Case 1
set(gcf,'Position',[0.05,0.5,0.95,0.4])  %2-D Waves - Case 3 (modified)
%set(gcf,'Position',[0.1,0.5,0.8,0.4])   %Original


%INITIALIZE PARTICLE SIZE, PLOT AXES AND
%POSITION OF TEXT
ParticleSize=12;
xmin=-0.05*vlx;
xmax=1.05*vlx;
zmin=-0.01;
zmax=1.5*vlz;

% LOAD AND PLOT IPART
PART=load('IPART');
if(nbf<nb)
    %   PLOT FIXED AND MOVING BOUNDARY PARTICLES AND WATER PARTICLES
    figure(XF)
    plot(PART(1:nbf,1),PART(1:nbf,2),'r.', 'MarkerSize',ParticleSize);     
    hold on;
    plot(PART(nbf+1:nb,1),PART(nbf+1:nb,2),'k.', 'MarkerSize',ParticleSize); 
    plot(PART(nb+1:np,1),PART(nb+1:np,2),'b.', 'MarkerSize',ParticleSize); 
else
    %   PLOT FIXED AND MOVING BOUNDARY PARTICLES AND WATER PARTICLES
    figure(XF)
    plot(PART(1:nb,1),PART(1:nb,2),'r.', 'MarkerSize',ParticleSize);     
    hold on;
    plot(PART(nb+1:np,1),PART(nb+1:np,2),'b.', 'MarkerSize',ParticleSize); 
end
%   SET AXES, LABELS AND TITLE
axis([xmin xmax zmin zmax])
xlabel('X(m)')
ylabel('Z(m)')
title('Initial Particle Configuration')