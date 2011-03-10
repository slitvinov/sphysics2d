.SUFFIXES: .f .obj .exe

OPTIONS=

OBJFILES=

.f.obj:
         ftn95 $< $(OPTIONS)

.obj.exe:
         slink $(OBJFILES) -FILE:$@

