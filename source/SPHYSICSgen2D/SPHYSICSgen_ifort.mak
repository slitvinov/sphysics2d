FC=ifort -O3
srcdir = .
idir=../../execs
bakdir=../../execs.bak
objects=SPHYSICSgen_2D.o \
#
SPHYSICSgen: $(objects)
	$(FC) -o SPHYSICSgen_2D $(srcdir)/$(objects)
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
	if [ -f $(idir)/SPHYSICSgen_2D ]; then \
	mv $(idir)/SPHYSICSgen_2D $(bakdir)/; \
	echo Old SPHYSICSgen_2D moved to execs.bak from execs; \
	fi
#
	mv SPHYSICSgen_2D $(idir)
	echo New SPHYSICSgen_2D moved to execs
#
clean :
	rm *.o
