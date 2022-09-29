.SUFFIXES: .f90 .o .mod

#default
PROGRAM     = bgk_thermal 
FC          = gfortran -cpp
CC          = gcc
COPT        = -O3
!!COPT        = -g -fbounds-check -fbacktrace
!!POPT        = 
!!POPT        = -fopenmp
FIX         =

ifdef INTEL
FC=ifx -fpp
FC=ifort -fpp
CC=icc 
COPT = -O3 -mtune=skylake-avx512
endif

ifdef NV
FC=nvfortran  -cpp
CC=nvcc 
COPT = -O3 -Minfo=all
endif


OBJS = \
BC_f.o \
BC_g.o \
phase_field.o \
forces.o \
regular.o \
collision.o \
equilibrium.o \
streaming.o \
moments.o \
moments2.o \
momentsPre.o \
dump_in.o \
dump_out.o \
probe.o \
prof_x.o \
prof_y.o \
setup.o \
initialization.o \
out2d.o \
mass_check.o \
media.o \
diag.o \
get_mem.o \
main.o \
##output_paraview_structured.o\

all	: compile link
compile : shared.o get_mem_mod.o $(OBJS)

$(OBJS) : 
.c.o  :
	$(CC) $(POPT) -c $(FIX) $<
.f90.o  :
	$(FC) $(POPT) -c $(COPT) $(FIX) $<
.mod.o  :
	$(FC) $(POPT) -c $(COPT) $(FIX) $<

link    :
	$(FC) $(POPT) $(FIX) -o RUN/$(PROGRAM) shared.o $(OBJS)
#	mv *.o RUN/
rmobj   :
	rm -f *.o *.mod *.vts *.vtu *.vtm *.gz dumpFile.dat log.dat parameters.dat CpCf.dat info_* bgk_thermal  

clean   :
	rm -f *.o *.mod
remake  : clean rmobj  all
