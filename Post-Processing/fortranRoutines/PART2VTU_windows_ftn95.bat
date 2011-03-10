set UDIRX= %CD%

cd ..\..\Post-Processing\fortranRoutines


mk32 -f PART2VTU_ftn95.mak clean
mk32 -f PART2VTU_ftn95.mak

copy PART2VTU_2D.exe ..\..\execs\

cd %UDIRX%


del Paraview*


md ParaviewFiles\VTU


copy ..\..\Post-Processing\fortranRoutines\PART2VTU_2D.exe PART2VTU_2D.exe

PART2VTU_2D.exe 











