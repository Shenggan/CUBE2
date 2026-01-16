#OPTIONS+=-CB
OPTIONS+=-DPID
OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o))) Green.o

all: main.x
	@echo "done"
main.x: $(OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(OBJFILE) -o $@ $(FFTFLAG)

$(OBJFILE): variables.o
parameters.o: Makefile basic_functions.f08
variables.o: parameters.o
pencil_fft.o: parameters.o
particle_mesh.o: variables.o pencil_fft.o cubefft.o
initialize.o: variables.o Green.o cubefft.o pencil_fft.o z_checkpoint.txt z_halofind.txt
finalize.o: variables.o cubefft.o pencil_fft.o
main.o: $(OBJFILE)
#$(OBJFILE): variables.o

parameters.o: parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
Green.o: ./Green/Green.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

clean:
	rm -f *.mod *.o *.out *.err *.x *~
