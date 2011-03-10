FC=ifort -O3
srcdir = .
idir=../../execs
bakdir=../../execs.bak
objects=PART2VTU_2D.o \
#
PART2VTU: $(objects)
	$(FC) -o PART2VTU_2D $(srcdir)/$(objects)
#
	if [ -d $(bakdir) ]; then \
	echo execs.bak Directory Exists; else \
	mkdir $(bakdir); \
	echo execs.bak Directory Created; \
	fi
#
	if [ -d $(idir) ]; then \
	echo execs Directory Exists; else \
	mkdir $(idir); \
	echo execs Directory Created; \
	fi

#
	if [ -f $(idir)/PART2VTU_2D ]; then \
	mv $(idir)/PART2VTU_2D $(bakdir)/; \
	echo Old PART2VTU_2D moved to execs.bak from execs; \
	fi
#
	mv PART2VTU_2D $(idir)
	echo New PART2VTU_2D moved to execs
#
clean :
	rm *.o
