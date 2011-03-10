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

%LOAD AND READ TIME,DT FROM FILE DT. THE FIRST HEADERLINE IN THE FILE IS
%SKIPPED
irun=input('Is the run complete? (1=Yes,0=No)  ');
if(irun==1)
    %LOAD AND READ TIME,DT FROM FILE DT. THE FIRST HEADLINE IN THE FILE IS
    %SKIPPED
    [time,DT1,DT2,DT]=textread('DT','%f %f %f %f','headerlines',1);
    %THE NUMBER OF PARTFILES, Nframes, IS ONE LESS THAN THE LENGTH OF THE time VECTOR.
    %THIS IS BECAUSE THE LAST VALUE OF THE TIME VECTOR IS ASSOCIATED WITH EPART
    Nframes=length(time)-1;
else
    Nframes=input('Number of PART files to be visualized= ');
end

%CHECK FOR EXISTENCE OF ParaviewFiles/VTU DIRECTORY. THIS IS THE DIRECTORY WHERE
% THE .vtu FILRS WILL BE STORED. IF IT DOES NOT
%EXIST CREATE IT
dirname='ParaviewFiles/VTU';
dirstatus=exist(dirname,'dir')
if(dirstatus==0)
    mkdir(dirname)
end
    
%THIS LOOP DETERMINES THE NUMBER OF DESIRED SUBDIVISIONS OF THE PARTICLE DATA FOR COLORING. THEN A UNIQUE SCALAR VALUE TITLED Scalarplot
% IS ASSIGNED FOR EACH DIVISION. THE SCALAR VALUE FOR DIFFERENT DIVISIONS ARE DIFFERENT. THIS IS ARBITRARILY CHOSEN
% SO AS TO INCLUDE THE ABILITY TO VIEW THE PARTICLES IN DIFFERENT DIVISIONS WITH DIFFERENT COLORS IN PARAVIEW
if(nbf<nb)
    idiv=3;
    nbeg1=1
    nend1=nbf;
    nbeg2=nbf+1;
    nend2=nb;
    nbeg3=nb+1
    nend3=np;    
else
    idiv=2;
    nbeg1=1
    nend1=nb;
    nbeg2=nb+1;
    nend2=np;
end
    
%LOOP OVER FRAMES    
    for iframe=1:Nframes
        iframe
        
% READ IN THE PART FILE FOR EACH FRAME
        if(iframe<10)
            eval(['PART=load(''PART_000',int2str(iframe),''');']);
            fid=eval(['fopen(''',dirname,'/PART00' int2str(iframe) '.vtu' ''',''w'' );']);
        elseif((iframe>=10)&(iframe<100))
            eval(['PART=load(''PART_00',int2str(iframe),''');']);
            fid=eval(['fopen(''',dirname,'/PART0' int2str(iframe) '.vtu' ''',''w'' );']);
        elseif((iframe>=100)&(iframe<1000))
            eval(['PART=load(''PART_0',int2str(iframe),''');']);
            fid=eval(['fopen(''',dirname,'/PART' int2str(iframe) '.vtu' ''',''w'' );']);
        else((iframe>=1000)&(iframe<10000))
            eval(['PART=load(''PART_',int2str(iframe),''');']);
            fid=eval(['fopen(''',dirname,'/PART' int2str(iframe) '.vtu' ''',''w'' );']);
        end
        
% READ POSITION, VELOCITY, DENSITY, PRESSURE, MASS AND VORTICITY DATA FOR ALL PARTICLES        
       
            xp=PART(:,1);
            zp=PART(:,2);            
            up=PART(:,3);            
            wp=PART(:,4);
            rhop=PART(:,5);
            P=PART(:,6);
            mass=PART(:,7);
       
    
%   OUTPUT TO FILE IN VTU FORMAT
        fprintf(fid,'<?xml version="1.0"?>\r\n');
        fprintf(fid,'<VTKFile type= "UnstructuredGrid"  version= "0.1"  byte_order= "BigEndian">\r\n');   
        fprintf(fid,' <UnstructuredGrid>\r\n');
        fprintf(fid,'  <Piece NumberOfPoints="%d" NumberOfCells="%d">\r\n',np,np);
        
% WRITE IN PRESSURE DATA        
        fprintf(fid,'   <PointData Scalars="Pressure" Vectors="Velocity">\r\n');
        fprintf(fid,'    <DataArray type="Float32" Name="Pressures" format="ascii">\r\n');
        for ii=1:np
            fprintf(fid,'%f\t',P(ii));
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        
% WRITE DENSITY DATA        
        fprintf(fid,'    <DataArray type="Float32" Name="Density" format="ascii">\r\n');
        for ii=1:np
            fprintf(fid,'%f\t',rhop(ii));
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        

        
% THIS SECTION IS USED TO COLOR DIFFERENT PARTICLES BASED THE INPUT IDIV SPECIFIED ABOVE.
        fprintf(fid,'    <DataArray type="Float32" Name="Scalarplot" format="ascii">\r\n');
        for ii=1:idiv
            eval(['nbeg=nbeg' int2str(ii) ';'])
            eval(['nend=nend' int2str(ii) ';'])
            for jj=nbeg:nend
                fprintf(fid,'%f\t',ii);
            end
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
       
% WRITE VELOCITY DATA
        fprintf(fid,'    <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\r\n');
        for ii=1:np
            vel=[up(ii) 0 wp(ii)];
            fprintf(fid,'%f\t %f\t %f\t',vel);
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        fprintf(fid,'   </PointData>\r\n');
        
% WRITE PARTICLE POSITION DATA
        fprintf(fid,'   <Points>\r\n');
        fprintf(fid,'    <DataArray type="Float32" NumberOfComponents="3" format="ascii">\r\n'); 
        for ii=1:np
            pos=[xp(ii) 0 zp(ii)];
            fprintf(fid,'%f\t %f\t %f\t',pos);
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        fprintf(fid,'   </Points>\r\n');
        
% WRITE CELL DATA. CELL IS OF TYPE VERTEX.        
        fprintf(fid,'   <Cells>\r\n');
        fprintf(fid,'    <DataArray type="Int32" Name="connectivity" format="ascii">\r\n');
        for ii=1:np
            fprintf(fid,'%d\t',ii-1);
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        fprintf(fid,'\r\n');
        fprintf(fid,'    <DataArray type="Int32" Name="offsets" format="ascii">\r\n');
        for ii=1:np
            fprintf(fid,'%d\t',ii);
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        fprintf(fid,'\r\n');
        fprintf(fid,'    <DataArray type="Int32" Name="types" format="ascii">\r\n');
        for ii=1:np
            fprintf(fid,'%d\t',1);
        end
        fprintf(fid,'\r\n');
        fprintf(fid,'    </DataArray>\r\n');
        fprintf(fid,'   </Cells>\r\n');   
        fprintf(fid,'  </Piece>\r\n')
        fprintf(fid,' </UnstructuredGrid>\r\n')
        fprintf(fid,'</VTKFile>')
    fclose(fid);
end
PVDinput