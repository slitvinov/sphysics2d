%THIS PROGRAM WRITES OUT A FILE CALLED VTUinp.pvd. THIS FILE IS NECESSARY TO ANIMATE THE VTU FILES (OBTAINED BY RUNNING PART2VTU.m) IN
% PARAVIEW 3. 
%  THIS FILE VTUinp.pvd SHOULD BE PLACED IN THE SAME DIRECTORY AS THE VTU FILES. THIS .pvd FILE BASICALLY CONTAINS A LIST OF VTU
% FILES THAT ARE TO BE ANIMATED. LOADING THIS FILE IN PARAVIEW 3.0 IS SUFFICIENT TO PERFORM THE ANIMATION


%OPEN THE VTUinp.pvd FILE FOR WRITING 
fid=fopen('ParaviewFiles/VTU/VTUinp.pvd','w')

% WRITE THE STANDARD HEADER FOR THE FILE    
fprintf(fid,'<?xml version="1.0"?>\r\n');
fprintf(fid,' <VTKFile type="Collection" version="0.1">\r\n');
fprintf(fid,'  <Collection>\r\n');

% WRITE THE FILE NAME FOR EACH TIME STEP
for ii=1:Nframes
    if ii<10
        fprintf(fid,'<DataSet timestep="%d" group="" part="%d" file="PART00%d.vtu"/>\r\n',ii,0,ii);
    elseif ii>99
        fprintf(fid,'<DataSet timestep="%d" group="" part="%d" file="PART%d.vtu"/>\r\n',ii,0,ii);
    else
        fprintf(fid,'<DataSet timestep="%d" group="" part="%d" file="PART0%d.vtu"/>\r\n',ii,0,ii);
    end 
end
    
fprintf(fid,'   </Collection>\r\n');
fprintf(fid,'  </VTKFile>\r\n');

fclose(fid);