# fortran compiler
FC   = mpif90
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
OBJ = initmpi.o types.o
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
common_mpi.o: common_mpi.f90 param.o
	$(FC) $(FFLAGS) -cpp -c  $<
initmpi.o: initmpi.f90 common_mpi.o types.o
	$(FC) $(FFLAGS) -cpp -c  $<
types.o: types.f90
	$(FC) $(FFLAGS) -cpp -c  $<
