#OPTIONS+=-CB
OPTIONS+=-DPID
OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o))) Green.o pp_force_kernel_cuda.o

XFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-num-images=1
OFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-num-images=1
FFTFLAG=-I${MKLROOT}/include/fftw/ -qmkl -L${CUDA_HOME}/lib64 -lcuda -lcudart -lstdc++ -latomic
# -L${CUDA_HOME}/lib64 -lcuda -lcudart -lstdc++

all: main.x
	@echo "done"
main.x: $(OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(OBJFILE) -o $@ $(FFTFLAG)

parameters.o: Makefile basic_functions.f08
variables.o: parameters.o
pencil_fft.o: parameters.o

# Fortran module ordering.  Keep this list explicit: making every object
# depend on variables.o also makes variables.o and parameters.o depend on
# themselves through OBJFILE, which creates circular dependencies.
buffer_grid.o buffer_particle.o checkpoint.o cubefft.o drift.o finalize.o \
initialize.o kick.o main.o particle_initialization.o particle_mesh.o timestep.o: variables.o

particle_mesh.o: variables.o pencil_fft.o cubefft.o
initialize.o: variables.o Green.o cubefft.o pencil_fft.o z_checkpoint.txt z_halofind.txt
finalize.o: variables.o cubefft.o pencil_fft.o

parameters.o: parameters.f90
	$(FC) -c $(OFLAG) $(OPTIONS) $<
Green.o: ./Green/Green.f90
	$(FC) -c $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile
	$(FC) -c $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

pp_force_kernel.o: pp_force_kernel.c
	icx -O3 -xHost -qopenmp -c pp_force_kernel.c -o pp_force_kernel.o
# 	gcc -O3 -fopenmp -c pp_force_kernel.c -o pp_force_kernel.o

pp_force_kernel_cuda.o: pp_force_kernel.cu
	${CUDA_HOME}/bin/nvcc -O3 -arch=sm_89 -use_fast_math -Xcompiler -fopenmp -c pp_force_kernel.cu -o pp_force_kernel_cuda.o

cuda_parameters_generated.h: parameters.f90 tools/gen_cuda_parameters.py
	python3 tools/gen_cuda_parameters.py parameters.f90 cuda_parameters_generated.h

pp_force_kernel_cuda.o: cuda_parameters_generated.h

clean:
	rm -f *.mod *.o *.out *.err *.x *~
