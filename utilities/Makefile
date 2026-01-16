OPTIONS=-DFFTFINE
OPTIONS+=-Dsigma8
OPTIONS+=-DPID
#OPTIONS+=-CB

MODEFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODEFILE:.f90=.o)))

all: ic.x cicpower.x fof.x fof_merge.x #write_gadget_format.x dsp.x
	@echo "made:" $^
ic.x: ic.o parameters.o pencil_fft.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
cicpower.x: cicpower.o parameters.o pencil_fft.o powerspectrum.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
dsp.x: parameters.o pencil_fft.o powerspectrum.o displacement.o
	$(FC) $(XFLAG) $(OPTIONS) $(FFTFLAG) $^ -o $@
write_gadget_format.x: write_gadget_format.o parameters.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
fof.x: parameters.o variables.o particle_initialization.o buffer_grid.o buffer_particle.o fof.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
fof_merge.x: parameters.o variables.o fof_merge.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@
multifields.x: parameters.o multifields.o
	$(FC) $(XFLAG) $(OPTIONS) $^ -o $@

parameters.o: Makefile ../basic_functions.f08
variables.o: parameters.o
particle_initialization.o: variables.o
buffer_grid.o: variables.o
buffer_particle.o: variables.o

parameters.o: ../parameters.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
variables.o: ../variables.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
particle_initialization.o: ../particle_initialization.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_grid.o: ../buffer_grid.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
buffer_particle.o: ../buffer_particle.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile parameters.o pencil_fft.o powerspectrum.o
	$(FC) $(OFLAG) $(OPTIONS) $(FFTFLAG) $< -o $@
pencil_fft.o: ../pencil_fft.f90 parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $< $(FFTFLAG)
*.o: Makefile

clean:
	rm -f *.mod *.o *.out *.err *.x *~
