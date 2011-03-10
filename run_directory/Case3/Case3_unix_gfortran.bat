UDIRX=`pwd`
cd ../../source/SPHYSICSgen2D/
make -f SPHYSICSgen_gfortran.mak clean
make -f SPHYSICSgen_gfortran.mak
if [ $? -eq 0 ]; then
  SPHYSICSGENcompilationDone="yes"
  echo 'SPHYSICSGENcompilationDone = ' $SPHYSICSGENcompilationDone
  cd $UDIRX
  ../../execs/SPHYSICSgen_2D < Case3.txt > "Case3.out"
  if [ $? -eq 0 ]; then
    cp SPHYSICS.mak ../../source/SPHYSICS2D/
    cd ../../source/SPHYSICS2D
    pwd
    make -f SPHYSICS_gfortran.mak clean
    make -f SPHYSICS_gfortran.mak 
    if [ $? -eq 0 ]; then
      SPHYSICScompilationDone="yes"
      echo 'SPHYSICScompilationDone = ' $SPHYSICScompilationDone
      rm SPHYSICS.mak
      cd $UDIRX
      pwd
      cp ../../source/SPHYSICS2D/SPHYSICS_2D SPHYSICS_2D
      ./SPHYSICS_2D
    else
      cd $UDIRX
      echo ' '
      echo 'SPHYSICS_2D compilation failed'
      echo 'Make sure correct compiler is selected in Case file'
    fi
  else
    echo ' '
    echo 'SPHYSICSgen_2D run failed'
    echo 'Make sure all parameters are correct in Case file'
  fi
else
  cd $UDIRX
  echo ' '
  echo 'SPHYSICSgen compilation failed'
fi
