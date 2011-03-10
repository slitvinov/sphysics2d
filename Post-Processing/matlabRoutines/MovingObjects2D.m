ll=0;
for ii=Nframes_ini:Nframes
    %READ IN EACH FRAME INTO PART1,PART2...PART'EFRAME'
    if ((ii>0)& (ii< 10)) 
        eval(['fn =','''PART_000', int2str(ii),'''',';'])
        PART=load(fn);
    elseif ((ii>=10) & (ii<100)) 
        eval(['fn =','''PART_00', int2str(ii),'''',';'])
        PART=load(fn);
    elseif ((ii>=100) & (ii<1000)) 
        eval(['fn =','''PART_0', int2str(ii),'''',';'])
        PART=load(fn);
    elseif ((ii>=1000) & (ii<10000)) 
        eval(['fn =','''PART_', int2str(ii),'''',';'])
        PART=load(fn);
    end
%   PLOT FIXED AND MOVING BOUNDARY PARTICLES AND WATER PARTICLES
    figure(XF)
    plot(PART(1:nbf,1),PART(1:nbf,2),'r.', 'MarkerSize',ParticleSize);     
    hold on;
    plot(PART(nbf+1:nb,1),PART(nbf+1:nb,2),'k.', 'MarkerSize',ParticleSize); 
    plot(PART(nb+1:np,1),PART(nb+1:np,2),'b.', 'MarkerSize',ParticleSize); 
%   SET AXES, LABELS AND TITLE
    axis([xmin xmax zmin zmax])
    xlabel('X(m)')
    ylabel('Z(m)')
    if(irun==1)
        title(['TIME=',num2str(time(ii)),'s : ','  FRAME=',num2str(ii)])
    else
        title(['  FRAME=',num2str(ii)])
    end
    hold off
    grid on
    h=get(gcf,'CurrentAxes');
    ll=ll+1;
%   GET FRAME FOR STORING A MOVIE    
    M1(:,ll)=getframe(gcf);
    hold off  
    clear PART 
    sframes = 0
    if(sframes==1)
    eval(['print -djpeg -r300 ', '''frame',num2str(ii),'.jpg''',';'])
    end
end
