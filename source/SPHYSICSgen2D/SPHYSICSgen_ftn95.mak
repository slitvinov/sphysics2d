## Uncomment/Comment the options as appropriate for debugging
OPTIONS=/OPTIMISE
#OPTIONS=/CHECK

OBJFILES=SPHYSICSgen_2D.obj

SPHYSICSgen_2D.exe: $(OBJFILES)

clean:
	del *.obj