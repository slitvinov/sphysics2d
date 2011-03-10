@ECHO OFF

set UDIRX= %CD%

cd ..\..\execs\

move PART2VTU_2D.exe ..\execs.bak

cd ..\Post-Processing\fortranRoutines\
del *.o

make -f PART2VTU_gfortran.mak

IF EXIST ../../execs/PART2VTU_2D.exe (
  ECHO.
  ECHO PART2VTU_2D compilation Done=yes
  ECHO.

  cd %UDIRX%
  
  IF EXIST ParaviewFiles\ (
    ..\..\execs\PART2VTU_2D.exe
    ) ELSE (
    md ParaviewFiles
    md ParaviewFiles\VTU
    ..\..\execs\PART2VTU_2D.exe
   )
  REM !!To use debugging in gfortran, edit PART2VTUgen_gfortran.mak

) ELSE (
  ECHO PART2VTU_2D compilation FAILED
  cd %UDIRX%
)









