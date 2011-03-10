FC=ifort

OPTIONS=/nologo
COPTIONS=/O3

OBJFILES=PART2VTU_2D.obj

.f.obj:
	$(FC) $(OPTIONS) $(COPTIONS) /c $<

PART2VTU_2D.exe: $(OBJFILES)
	xilink /OUT:$@ $(OPTIONS) $(OBJFILES)

clean :
	del *.obj *.mod PART2VTU_2D.exe
