# Compiling SNaC

For most systems, SNaC can be compiled from the root directory with `make`. The same CaNS-style target surface is available, so `make libs && make` is also accepted, although `make libs` is a no-op in SNaC because HYPRE and FFTW are expected to be provided by the user/system.

The `Makefile` in the root directory is used to compile the code. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), a few pre-defined profiles depending on the nature of the run (e.g. production vs debugging), and pre-processing options:

```shell
#
# compiler and compiling profile
#
FCOMP=GNU          # options: GNU, INTEL, INTEL_IFORT, NVIDIA, CRAY, FUJITSU
FFLAGS_OPT=1       # for production runs
FFLAGS_OPT_MAX=0   # for production runs (more aggressive optimization)
FFLAGS_DEBUG=0     # for debugging
FFLAGS_DEBUG_MAX=0 # for thorough debugging
#
# defines
#
SINGLE_PRECISION=0       # perform the whole calculation in single precision
TIMING=0                 # report wall-clock time per time step
IMPDIFF=0                # implicit discretization of the diffusion terms
BOUSSINESQ_BUOYANCY=0    # Boussinesq source term from the first scalar field
FFT_AXIS=3               # options: 0 (off), 1 (_FFT_X), 2 (_FFT_Y), 3 (_FFT_Z)
FFT_USE_SLABS=0          # use slab communication for the FFT path
FFT_USE_SLICED_PENCILS=0 # use sliced pencils for the FFT path
ONE_PRESS_CORR=0         # use one pressure correction when possible
OPENMP=0                 # OpenMP build
```

In this file, `FCOMP` can be one of `GNU`, `INTEL`, `INTEL_IFORT`, `NVIDIA`, `CRAY`, or `FUJITSU`; the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` options. Finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `configs/flags.mk`.

The following pre-processing options are available through `build.conf`:

 * `SINGLE_PRECISION`       : calculation will be carried out in single precision (the default precision is double)
 * `TIMING`                 : wall-clock time per time step is computed
 * `IMPDIFF`                : diffusion term of the N-S equations is integrated in time with an implicit discretization
 * `BOUSSINESQ_BUOYANCY`    : enables the Boussinesq source term from the first scalar field
 * `FFT_AXIS`               : `0` disables FFT acceleration; `1`, `2`, or `3` defines `_FFT_X`, `_FFT_Y`, or `_FFT_Z`, respectively
 * `FFT_USE_SLABS`          : enables the slab communication path for FFT acceleration
 * `FFT_USE_SLICED_PENCILS` : enables the sliced-pencil communication path for FFT acceleration
 * `ONE_PRESS_CORR`         : enables the one-pressure-correction path where supported
 * `OPENMP`                 : enables OpenMP compiler and linker flags

HYPRE is linked through `configs/libs.mk`. The default HYPRE library path is `$(HOME)/hypre/lib`; this can be overridden by editing `configs/libs.mk`, by setting `HYPRE_LIB_DIR`, or by passing, e.g., `make HYPRE_LIB_DIR=[PATH_TO_LIB]`. If HYPRE modules or headers are not found by the compiler wrapper, set `HYPRE_INC_DIR` as well.

FFTW is linked from the paths listed in `LD_LIBRARY_PATH`, with `pkg-config` used when available and `-lfftw3` as fallback. For single precision builds, `fftw3f` is linked as well. Additional linking and include options can be changed in `configs/libs.mk`, or by building with, e.g., `make LIBS+='-L[PATH_TO_LIB] -l[NAME_OF_LIB]' INCS+='-I[PATH_TO_INCLUDE]'`.

Typing `make` will compile the code and copy the executable `snac` to a `run/` folder. Typing `make run` will also copy the default input files `dns.nml` and `blocks.nml` under `src/`, plus any optional `src/grid/` folder, to the same `run/` folder.

`make clean` clears the SNaC build files, `make libsclean` is a no-op external-library cleanup target, and `make allclean` runs both. The default `build.conf` file and `*.mk` files are created from `configs/defaults` at the first compilation.
