FC=ifort

OPTIONS=/nologo
COPTIONS=/O3

OBJFILES=SPHYSICSgen_2D.obj

.f.obj:
	$(FC) $(OPTIONS) $(COPTIONS) /c $<

SPHYSICSgen_2D.exe: $(OBJFILES)
	xilink /OUT:$@ $(OPTIONS) $(OBJFILES)

clean :
	del *.obj *.mod SPHYSICSgen_2D.exe
