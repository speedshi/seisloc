# Compiler name
FC = gfortran
# MPI compiler name
MPIFC = mpif90
# OpenMP flag
OMP = -fopenmp
# Other compilation flags
FFLAGS = -g -O2 -ffree-line-length-none $(OMP)
# Option to set where module files are saved
MOD = -J

# Directory to executables
B = ../bin
# Object file directory
O = ../obj

# Module object files
MODS = \
    $(O)/paramod.o

# Other object files
OBJS = \
    $(O)/calcharfun.o \
    $(O)/combnum.o \
    $(O)/corrcoefm.o \
    $(O)/envhlbfft.o \
    $(O)/eventidf_3ds.o \
    $(O)/eventidf_stm.o \
    $(O)/evprocesst.o \
    $(O)/stkcorrcoef.o

# Programs
MPI_PROGS = \
    $(B)/migrationloc_MPI

OPENMP_PROGS = \
    $(B)/migrationloc_OpenMP \
    $(B)/migrationloc_OpenMP_utsp

SERIAL_PROGS = \
    $(B)/migrationloc_serial \

PROGS = $(MPI_PROGS) $(OPENMP_PROGS) $(SERIAL_PROGS)

# Default compilation rule: build all programs
default: $(PROGS)

mpi: $(MPI_PROGS)
omp: $(OPENMP_PROGS)
serial: $(SERIAL_PROGS)

## Compilation rules
# Modules must be compiled first in case object files depend on them
$(OBJS): $(MODS)

# Normal source files
$(O)/%.o: %.F90
	$(FC) $(FFLAGS) -c -I $(O) $(MOD) $(O) -o $@ $*.F90

# MPI programs
$(B)/%_MPI: %_MPI.F90 $(MODS) $(OBJS)
	$(MPIFC) $(FFLAGS) -I $(O) -o $@ $^

# Other programs
$(B)/%: %.F90 $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -I $(O) -o $@ $^

# Housekeeping rules
.PHONY: clean cleanall

clean:
	/bin/rm -f $(O)/*.o $(O)/*.mod

cleanall: clean
	/bin/rm -f $(PROGS)
