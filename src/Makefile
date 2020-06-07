# fortran compiler
FC   = mpifort
# flags
BIG :=
DBG := -O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all# -finit-real=snan -ffpe-trap=invalid 
OMP :=
OPT := #-O3
OTH := 
PROF := #-pg
FFLAGS := $(OPT) $(BIG) $(DBG) $(PROF) $(OMP) $(OTH)
TARGET = snac
#
OBJ = bound.o debug.o initflow.o initgrid.o chkdiv.o chkdt.o common_mpi.o correc.o fillps.o initmpi.o load.o main.o mom.o output.o param.o pressure_update.o rk.o scal.o sanity.o solver.o types.o
#
HYPRE_LIB_DIR=/home/pedro/software/hypre-2.11.2/lib
LIBS = -L$(HYPRE_LIB_DIR) -lHYPRE
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
	@printf "\nExecutable file $(TARGET) and input file dns.in copied to run folder $(RUNDIR)\n"
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
fillps.o: fillps.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initflow.o: initflow.f90 common_mpi.o param.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initgrid.o: initgrid.f90 param.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initmpi.o: initmpi.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
load.o: load.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
main.o: main.f90 bound.o debug.o initflow.o initgrid.o chkdiv.o chkdt.o common_mpi.o correc.o fillps.o initmpi.o load.o mom.o output.o param.o pressure_update.o rk.o scal.o sanity.o solver.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
mom.o: mom.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
output.o: output.f90 common_mpi.o load.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
param.o: param.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
pressure_update.o: pressure_update.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
rk.o: rk.f90 mom.o scal.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
scal.o: scal.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
sanity.o: sanity.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
solver.o: solver.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
types.o: types.f90
	$(FC) $(FFLAGS) -cpp -c  $<
