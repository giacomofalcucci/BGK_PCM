.SUFFIXES: .f90 .o .mod

PROGRAM     = bgk_thermal 
FC          = gfortran
CC          = gcc
COPT        = -O3
!!COPT        = -g -fbounds-check -fbacktrace
!!POPT        = 
!!POPT        = -fopenmp
FIX         =

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
dump_in.o \
dump_out.o \
setup.o \
initialization.o \
out2d.o \
mass_check.o \
media.o \
diag.o \
main.o \
##output_paraview_structured.o\

all	: compile link
compile : shared.o $(OBJS)

$(OBJS) : 
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
