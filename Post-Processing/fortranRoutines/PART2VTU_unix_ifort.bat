UDIRX=`pwd`
cd ../../Post-Processing/fortranRoutines/
make -f PART2VTU_ifort.mak clean
make -f PART2VTU_ifort.mak
cd $UDIRX
	if [ -d ./ParaviewFiles/VTU ]; then \
	echo ParaviewFiles Directory Exists; else \
	mkdir ParaviewFiles; \
        mkdir ParaviewFiles/VTU; \
	echo ParaviewFiles Directory Created; \
	fi
../../execs/PART2VTU_2D
