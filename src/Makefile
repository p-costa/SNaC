# fortran compiler
FC   = mpifort
# flags
BIG :=
DBG := -g -fbacktrace -Wall -O0 -Wextra -pedantic -fcheck=all
OMP :=
OPT := -O3
OTH := 
PROF := #-pg
FFLAGS := $(OPT) $(BIG) $(DBG) $(PROF) $(OMP) $(OTH)
TARGET = cans
#
OBJ = bound.o initgrid.o chkdiv.o chkdt.o common_mpi.o correc.o fillps.o initmpi.o load.o mom.o param.o scal.o types.o
#
LIBS =
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
bound.o: bound.f90 types.o
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
initgrid.o: initgrid.f90 param.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
initmpi.o: initmpi.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
load.o: load.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
mom.o: mom.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
param.o: param.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
scal.o: scal.f90 types.o
	$(FC) $(FFLAGS) -cpp -c  $<
types.o: types.f90
	$(FC) $(FFLAGS) -cpp -c  $<
