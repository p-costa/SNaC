## Synopsis

**SNaC** is [**CaNS**](https://github.com/p-costa/CaNS) spelled backwards, and is a multi-block code for massively parallel direct numerical simulations (DNS) of fluid flows. *SNaC* aims at combining the versatility of a multi-block DNS solver, with the FFT-based acceleration used in CaNS.

The solver is able to simulate the flow in any three-dimensional multi-block structured Cartesian grid. However, if the geometry has one homogeneous, 'extruded' direction with constant grid spacing, SNaC can use a very fast solver that exploits FFTs in this direction. This is SNaC's *warp drive* :rocket:, as it yields a huge speedup in wall-clock time per time step.

**Reference**

P. Costa. *A FFT-accelerated multi-block finite-difference solver for massively parallel simulations of incompressible flows*.
Comput. Phys. Commun. 271 : 108194 (2022) [[DOI:10.1016/j.cpc.2021.108194]](doi.org/10.1016/j.cpc.2021.108194) [[arXiv:2106.03583]](https://arxiv.org/pdf/2106.03583.pdf).

## News
[08/07/2022] The input files describing the block geometry (under `geo/block.???`) have been simplified. Now, instead of prescribing the lower and upper extents of each block `lo(:)` and `hi(:)`, the number of grid points `ng(:)` is prescribed. This change makes it much easier to refine the grid, since one does not need to correct extents of neighboring blocks. See the updated [`src/INFO_INPUT.md`](src/INFO_INPUT.md) for more details.

## Features

Some features are:

 * Multi-block, three-dimensional parallelization
 * Hybrid MPI/OpenMP parallelization
 * FFT-based synthesis of the Poissonn equation along one direction
 * HYPRE library used to solve Poisson/Helmholtz equations
 * Parallel I/O using MPI I/O
 * A different canonical flow can be simulated just by changing the input files

## Motivation

*SNaC* is meant to serve as a multi-block DNS code for fast, massively-parallel simulations of single-phase flows, and as a solid base solver on top of which more complex phenomena can be implemented, such as numerical methods for multiphase flows.

## Method

The fluid flow is solved with a standard second-order finite-difference/-volume pressure correction scheme. Time is advanced with a three-step low storage Runge-Kutta scheme. Optionally, for increased stability at low Reynolds numbers, at the price of higher computational demand, the diffusion term can be treated implicitly.

## Usage

### Input files

The input files `dns.in` sets the physical and computational parameters, while the block files `geo/block.???` setup block-specific parameters. In the `examples/` folder are examples of input files for several canonical flows. See [`src/INFO_INPUT.md`](src/INFO_INPUT.md) for a detailed description of the input files.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2-, and 3-dimensional output files, respectively. *The code should be recompiled after editing out?d.h90 files*.

### Build

The code should be compiled in `src/`. The prerequisites are the following:

 * MPI
 * [*HYPRE*](https://github.com/hypre-space/hypre)
 * OpenMP (optional)
 * *FFTW* (optional, in case FFT acceleration is used)

The Makefile in `src/` should be modified in agreement to the installation paths of each library. Also, the following preprocessor options are available:

 * `-D_TIMING`           : wall-clock time per time step is computed
 * `-D_IMPDIFF`          : diffusion term of the N-S equations is integrated in time with an implicit discretization (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `-D_SINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double)
 * `-D_FFT_?`, with `?` being `X`, `Y` or `Z`: will use FFTs to solve the Poisson equation in the direction in question.

Typing `make run` will compile the code and copy the executable `snac` and input file `dns.in` to a `run/` folder.

### Running the code

Run the executable with `mpirun` with a number of tasks and shared threads complying to what has been set in the input file `dns.in` (or in the `geo/block.???` files in case of multi-block). Data will be written by default in a folder named `data/`, which must be located where the executable is run.

### Visualizing field data

See [`src/INFO_VISU.md`](src/INFO_VISU.md).

## Notes

I appreciate any feedback that can improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.

Please read the `LICENSE` file.

## Contributors

Pedro Costa (p.simoes.costa@gmail.com)
