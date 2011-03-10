@ECHO OFF
del *.exe
del *.mak

set UDIRX= %CD%

cd ..\..\execs\

move SPHYSICSgen_2D.exe ..\execs.bak

cd ..\source\SPHYSICSgen2D
del *.exe

nmake -f SPHYSICSgen_win_ifort.mak clean
nmake -f SPHYSICSgen_win_ifort.mak

IF EXIST SPHYSICSgen_2D.exe (
  ECHO.
  ECHO SPHYSICSGEN compilation Done=yes
  ECHO.
  copy SPHYSICSgen_2D.exe ..\..\execs\SPHYSICSgen_2D.exe

  cd %UDIRX%

  copy ..\..\execs\SPHYSICSgen_2D.exe SPHYSICSgen_2D.exe

  SPHYSICSgen_2D.exe <Case1.txt > Case1.out

  copy SPHYSICS.mak ..\..\source\SPHYSICS2D\SPHYSICS.mak

  cd ..\..\execs\

  del *.obj

  move SPHYSICS_2D.exe ..\execs.bak

  cd ..\source\SPHYSICS2D
  del *.exe

  nmake -f SPHYSICS.mak clean
  nmake -f SPHYSICS.mak

  IF EXIST SPHYSICS_2D.exe (
    ECHO.
    ECHO SPHYSICScompilationDone = yes
    ECHO.
    copy SPHYSICS_2D.exe ..\..\execs\SPHYSICS_2D.exe

    cd %UDIRX%

    copy ..\..\execs\SPHYSICS_2D.exe SPHYSICS_2D.exe 

    SPHYSICS_2D.exe

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









