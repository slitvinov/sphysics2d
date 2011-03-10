@ECHO OFF
del *.exe
del *.mak

set UDIRX= %CD%

cd ..\..\execs\

move SPHYSICSgen_2D.exe ..\execs.bak

cd ..\source\SPHYSICSgen2D
del *.exe

mk32 -f SPHYSICSgen_ftn95.mak clean
mk32 -f SPHYSICSgen_ftn95.mak

IF EXIST SPHYSICSgen_2D.exe (
  ECHO.
  ECHO SPHYSICSGEN compilation Done=yes
  ECHO.
  copy SPHYSICSgen_2D.exe ..\..\execs\SPHYSICSgen_2D.exe

  cd %UDIRX%

  copy ..\..\execs\SPHYSICSgen_2D.exe SPHYSICSgen_2D.exe

  SPHYSICSgen_2D.exe <Case1.txt > Case1.out
  REM -- To use debugging in ftn95, use following line instead (remove REM, and comment out above line with REM)
  REM sdbg SPHYSICSgen_2D.exe <Case1.txt > Case1.out
  REM -- NOTE: you must compile SPHYSICSgen_2D.f with the /CHECK compiling option in SPHYSICSgen_ftn95.mak

  copy SPHYSICS.mak ..\..\source\SPHYSICS2D\SPHYSICS.mak

  cd ..\..\execs\

  del *.obj

  move SPHYSICS_2D.exe ..\execs.bak

  cd ..\source\SPHYSICS2D
  del *.exe

  mk32 -f SPHYSICS.mak clean
  mk32 -f SPHYSICS.mak

  IF EXIST SPHYSICS_2D.exe (
    ECHO.
    ECHO SPHYSICScompilationDone = yes
    ECHO.
    copy SPHYSICS_2D.exe ..\..\execs\SPHYSICS_2D.exe

    cd %UDIRX%

    copy ..\..\execs\SPHYSICS_2D.exe SPHYSICS_2D.exe 

    SPHYSICS_2D.exe
    REM -- To use debugging in ftn95, use following line instead (remove REM, and comment out above line with REM)
    REM sdbg SPHYSICS_2D.exe
    REM -- NOTE: you must compile SPHYSICS_2D.exe with the /CHECK compiling option,  
    REM -- EDIT the 'tocompile_ftn95' subroutine of SPHYSICSgen_2D.f

  ) ELSE (
    ECHO.
    ECHO SPHYSICScompilation FAILED
    ECHO Check you have specified the correct compiler in Case file
    ECHO.

    cd %UDIRX%
  )
) ELSE (
  ECHO SPHYSICSGEN compilation FAILED
  cd %UDIRX%
)









