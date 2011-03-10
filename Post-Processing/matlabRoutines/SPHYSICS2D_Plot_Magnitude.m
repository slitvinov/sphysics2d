clear;

load matlabin -ASCII
np=matlabin(1,1)
bxmax=matlabin(2,1)
bzmax=matlabin(4,1)
out=matlabin(5,1);
nb=matlabin(6,1);

load INDAT -ASCII
tmax=INDAT(31,1)
out=INDAT(32,1)

load IPART -ASCII
xmin=min(IPART(1:nb,1));
xmax=max(IPART(1:nb,1));
zmin=min(IPART(1:nb,2));
zmax=max(IPART(1:nb,2));

ini=input('Initial image ??');
ien=input('Last image ??');
kind=input('Property to plot: density (5), pressure (6), mass (7), vorticity (8) ??');
limit=input('Limits: manuals (1) or automatics (2) ??');
step_plot=input('Step for valticks ?? ');
if limit==1
lim_inf=input('   Limit inf ??');
lim_sup=input('   Limit sup ??'); 
end
ini_p=input('First particle to plot ??');
if ini_p==0
    ini_p=1
end
i_save=input('Save the images (y=1) ??');

for i=ini:ien
    %figure(1)
    axis equal;
    clf;
    colormap=jet;

    if i <10
       name=sprintf('PART_000%d',i)
    end
    if i>= 10 & i<=100
       name=sprintf('PART_00%d',i)
    end
    if i>= 100
       name=sprintf('PART_0%d',i)
    end

    eval([ '!copy ' name ' PART'])
  
    load PART  -ASCII;

    %istart1=1;
    %istop1=istart1+nb;
    %istop=np; 
    axis([xmin xmax zmin zmax]);  
    hold on;       
    VAR=PART(ini_p:np,kind);
        
    if limit==2
        clim_inf=min(VAR);
        clim_sup=max(VAR);
        lim_inf=clim_inf;
        lim_sup=clim_sup;       
    end 
      
    for j=ini_p:np %Color Scale
        variable=PART(j,kind);        
        PcolorIndex=round((variable-lim_inf)*64/(lim_sup-lim_inf));
     %    fprintf(1,'%d %f %d \n',j,vorticity,PcolorIndex); 
        if PcolorIndex > 64
            PcolorIndex = 64;
        end
        if PcolorIndex < 1
            PcolorIndex = 1;
        end
        Pcolor = colormap(PcolorIndex,:);
        newP = plot(PART(j,1),PART(j,2),'.','MarkerFaceColor', Pcolor);  
        hold on;
        set(newP,'Color',Pcolor); 
   end 
 
   hold on;
   cb = colorbar('horiz');
   %cb = colorbar;
   revalticks=lim_inf:step_plot:lim_sup;
   valticks=round((revalticks-lim_inf)*64/(lim_sup-lim_inf));
   set(cb,'Xtick',[valticks]);
   set(cb,'Xticklabel',[revalticks]);
     
   hold on;
   %time=out*i;
   %texto=sprintf('T= %4.2f s',time);
   %text(0.4,0.9,texto);
   xlabel(i);   
     
   hold off;

   if i_save==1
        eval(['print -djpeg ' name '.jpeg'])
   end

end
