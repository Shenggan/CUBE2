#OPTIONS+=-CB
OPTIONS+=-DPID
OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o))) Green.o pp_force_kernel.o

XFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-num-images=1
OFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-num-images=1
FFTFLAG=-I${MKLROOT}/include/fftw/ -qmkl

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
	$(FC) -c $(OFLAG) $(OPTIONS) $<
Green.o: ./Green/Green.f90
	$(FC) -c $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile
	$(FC) -c $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

pp_force_kernel.o: pp_force_kernel.c
	icx -O3 -xHost -qopenmp -c pp_force_kernel.c -o pp_force_kernel.o
# 	gcc -O3 -fopenmp -c pp_force_kernel.c -o pp_force_kernel.o

clean:
	rm -f *.mod *.o *.out *.err *.x *~
