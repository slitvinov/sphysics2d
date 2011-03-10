@ECHO OFF
del *.exe
del *.mak

set UDIRX= %CD%

cd ..\..\execs\

move SPHYSICSgen_2D.exe ..\execs.bak

cd ..\source\SPHYSICSgen2D
del *.exe

make -f SPHYSICSgen_gfortran.mak clean
make -f SPHYSICSgen_gfortran.mak

IF EXIST ../../execs/SPHYSICSgen_2D.exe (
  ECHO.
  ECHO SPHYSICSGEN compilation Done=yes
  ECHO.
  copy SPHYSICSgen_2D.exe ..\..\execs\SPHYSICSgen_2D.exe

  cd %UDIRX%

  copy ..\..\execs\SPHYSICSgen_2D.exe SPHYSICSgen_2D.exe

  SPHYSICSgen_2D.exe <Case4.txt > Case4.out
  REM !!To use debugging in gfortran, edit SPHYSICSgen_gfortran.mak

  copy SPHYSICS.mak ..\..\source\SPHYSICS2D\SPHYSICS.mak

  cd ..\..\execs\

  move SPHYSICS_2D.exe ..\execs.bak\

  cd ..\source\SPHYSICS2D
  del *.exe

  make -f SPHYSICS.mak clean
  make -f SPHYSICS.mak

  IF EXIST ../../execs/SPHYSICS_2D.exe (
    ECHO.
    ECHO SPHYSICScompilationDone = yes
    ECHO.
    copy SPHYSICS_2D.exe ..\..\execs\SPHYSICS_2D.exe

    cd %UDIRX%

    copy ..\..\execs\SPHYSICS_2D.exe SPHYSICS_2D.exe 

    SPHYSICS_2D.exe
    REM -- To use debugging in gfortran, edit subroutine 'tocompile_gfortran' in SPHYSICSgen_2D.f

  ) ELSE (
    ECHO.
    ECHO SPHYSICS compilation FAILED
    ECHO Check you have specified the correct compiler in Case file
    ECHO.

    cd %UDIRX%
  )
) ELSE (
  ECHO SPHYSICSGEN compilation FAILED
  cd %UDIRX%
)









