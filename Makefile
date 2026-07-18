#OPTIONS+=-CB
OPTIONS+=-DPID
OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
CPU_OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.o))) Green.o
CUDA_OBJFILE:=$(addprefix ,$(notdir $(MODFILE:.f90=.cuda.o))) Green.cuda.o pp_force_kernel_cuda.o pm3_fft_cuda.o
CUDA_MODDIR:=.cuda-mod

XFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-num-images=1
OFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large -coarray-num-images=1
FFTFLAG=-I${MKLROOT}/include/fftw/ -qmkl -latomic
CUDA_LDFLAGS=-L${CUDA_HOME}/lib64 -lcuda -lcudart -lcufft -lstdc++
CUDA_FFLAGS=$(OFLAG) $(OPTIONS) -DCUDA -module $(CUDA_MODDIR) -I$(CUDA_MODDIR)

all: main.x
	@echo "done"
main.x: $(CPU_OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(CPU_OBJFILE) -o $@ $(FFTFLAG)

.PHONY: all cuda clean

# Build the CUDA PP kernel and compile the Fortran sources with the CUDA macro.
# Separate object and module names let CPU and CUDA builds coexist safely.
cuda: main.cuda.x

main.cuda.x: $(CUDA_OBJFILE)
	@echo "Link CUDA files:"
	$(FC) $(XFLAG) $(OPTIONS) -DCUDA $(CUDA_OBJFILE) -o $@ $(FFTFLAG) $(CUDA_LDFLAGS)

parameters.o: Makefile basic_functions.f08
variables.o: parameters.o
pencil_fft.o: parameters.o

parameters.cuda.o: Makefile basic_functions.f08
variables.cuda.o: parameters.cuda.o
pencil_fft.cuda.o: parameters.cuda.o

# Fortran module ordering.  Keep this list explicit: making every object
# depend on variables.o also makes variables.o and parameters.o depend on
# themselves through OBJFILE, which creates circular dependencies.
buffer_grid.o buffer_particle.o checkpoint.o cubefft.o drift.o finalize.o \
initialize.o kick.o main.o particle_initialization.o particle_mesh.o timestep.o: variables.o

buffer_grid.cuda.o buffer_particle.cuda.o checkpoint.cuda.o cubefft.cuda.o drift.cuda.o finalize.cuda.o \
initialize.cuda.o kick.cuda.o main.cuda.o particle_initialization.cuda.o particle_mesh.cuda.o timestep.cuda.o: variables.cuda.o

particle_mesh.o: variables.o pencil_fft.o cubefft.o
initialize.o: variables.o Green.o cubefft.o pencil_fft.o z_checkpoint.txt z_halofind.txt
finalize.o: variables.o cubefft.o pencil_fft.o

particle_mesh.cuda.o: variables.cuda.o pencil_fft.cuda.o cubefft.cuda.o
initialize.cuda.o: variables.cuda.o Green.cuda.o cubefft.cuda.o pencil_fft.cuda.o z_checkpoint.txt z_halofind.txt
finalize.cuda.o: variables.cuda.o cubefft.cuda.o pencil_fft.cuda.o

parameters.o: parameters.f90
	$(FC) -c $(OFLAG) $(OPTIONS) $<
Green.o: ./Green/Green.f90
	$(FC) -c $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile
	$(FC) -c $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

$(CUDA_MODDIR):
	mkdir -p $@

Green.cuda.o: ./Green/Green.f90 | $(CUDA_MODDIR)
	$(FC) -c $(CUDA_FFLAGS) $< -o $@
%.cuda.o: %.f90 Makefile | $(CUDA_MODDIR)
	$(FC) -c $(CUDA_FFLAGS) $< -o $@ $(FFTFLAG)

pp_force_kernel_cuda.o: pp_force_kernel.cu
	${CUDA_HOME}/bin/nvcc -O3 -arch=sm_89 -use_fast_math -DCUDA -Xcompiler -fopenmp -c pp_force_kernel.cu -o pp_force_kernel_cuda.o

pm3_fft_cuda.o: pm3_fft_cuda.cu cuda_parameters_generated.h
	${CUDA_HOME}/bin/nvcc -O3 -arch=sm_89 -c pm3_fft_cuda.cu -o pm3_fft_cuda.o

cuda_parameters_generated.h: parameters.f90 tools/gen_cuda_parameters.py
	python3 tools/gen_cuda_parameters.py parameters.f90 cuda_parameters_generated.h

pp_force_kernel_cuda.o: cuda_parameters_generated.h

clean:
	rm -f *.mod *.o *.out *.err *.x *~
	rm -rf $(CUDA_MODDIR)
