# fortran compiler
FC   = mpifort
# flags
BIG :=
DBG := #-O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
OMP := #-fopenmp
OPT := -O3
OTH := -D_FFT_  -D_NON_NEWTONIAN -D_IMPDIFF
PROF := #-pg
FFLAGS := $(OPT) $(BIG) $(DBG) $(PROF) $(OMP) $(OTH)
TARGET = snac
#
OBJ = bound.o debug.o initflow.o initgrid.o chkdiv.o chkdt.o common_mpi.o correc.o fftw.o fft.o fillps.o initmpi.o load.o main.o mom.o non_newtonian.o output.o param.o updt_pressure.o rk.o scal.o sanity.o solver.o types.o
#
HYPRE_LIB_DIR=$(HOME)/hypre/lib
LIBS = -L$(HYPRE_LIB_DIR) -lHYPRE -lfftw3 -lfftw3_threads
#
all: $(TARGET)
#
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o $(TARGET)
# src, run and data directories
SRCDIR := $(pwd)
RUNDIR := $(SRCDIR)../run
DATDIR := $(RUNDIR)/data
#
run: $(TARGET)
	@mkdir -p $(RUNDIR) $(DATDIR)
	@cp $(TARGET) $(RUNDIR)
	@cp dns.in $(RUNDIR)
	@cp -r geo $(RUNDIR)
	@printf "\nExecutable file $(TARGET), input file dns.in and geo folder copied to run folder $(RUNDIR)\n"
#
.PHONY: clean
clean:
	rm -rf *.o *.mod *dSYM $(TARGET)
# dependencies
bound.o: bound.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
debug.o: debug.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
chkdiv.o: chkdiv.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
chkdt.o: chkdt.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
common_mpi.o: common_mpi.f90
	$(FC) $(FFLAGS) -cpp -c  $<
correc.o: correc.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
fft.o: fft.f90 common_mpi.o fftw.o param.o
	$(FC) $(FFLAGS) -cpp -c  $<
fftw.o: fftw.f90
	$(FC) $(FFLAGS) -cpp -c  $<
fillps.o: fillps.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initflow.o: initflow.f90 common_mpi.o param.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initgrid.o: initgrid.f90 common_mpi.o param.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initmpi.o: initmpi.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
load.o: load.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
main.o: main.f90 bound.o debug.o initflow.o initgrid.o chkdiv.o chkdt.o common_mpi.o correc.o fft.o fillps.o initmpi.o load.o mom.o non_newtonian.o output.o param.o updt_pressure.o rk.o scal.o sanity.o solver.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
mom.o: mom.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
non_newtonian.o: non_newtonian.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
output.o: output.f90 common_mpi.o load.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
param.o: param.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
updt_pressure.o: updt_pressure.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
rk.o: rk.f90 non_newtonian.o mom.o scal.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
scal.o: scal.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
sanity.o: sanity.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
solver.o: solver.f90 common_mpi.o fft.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
types.o: types.f90
	$(FC) $(FFLAGS) -cpp -c  $<
