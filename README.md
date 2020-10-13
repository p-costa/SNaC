## Synopsis

**SNaC** is [**CaNS**](https://github.com/p-costa/CaNS) spelled backwards, and is a code for massively-parallel direct numerical simulations (DNS) of fluid flows. *SNaC* is an alias for a longer and somewhat malevolent name: *Slow CaNS*. But what this tool loses in terms of speed compared to CaNS, it gains in versatility.

Indeed, a **multi-block** implementation of *SNaC* has been recently developed in branch [`multi_block`](https://github.com/p-costa/SNaC/tree/multi_block). Moreover, if the geometry has one homogeneous direction with constant grid spacing, branch [`multi_block_fft`](https://github.com/p-costa/SNaC/tree/multi_block) implements a very fast solver that exploits FFTs in this direction. This is SNaC's *warp drive* :rocket:, as yields a huge speedup in wall-clock time per time step.

## News

A **multi-block** implementation of *SNaC* has been developed in branch [`multi_block`](https://github.com/p-costa/SNaC/tree/multi_block) (and more recently [`multi_block_fft`](https://github.com/p-costa/SNaC/tree/multi_block)).

## Features

Some features are:

 * Block-structured three-dimensional parallelization
 * Hybrid MPI/OpenMP parallelization
 * HYPRE library used to solve Poisson/Helmholtz equations
 * Parallel I/O using MPI I/O 
 * A different canonical flow can be simulated just by changing the input files
 * Simple implementation with no dependencies on external libraries other than HYPRE and MPI

## Motivation

This code for canonical fluid flows will eventually serve as base for a multi-block DNS code for multiphase flows with adaptive mesh refinement.

## Method

The fluid flow is solved with a second-order finite-volume pressure correction scheme, discretized in a MAC grid arrangement. Time is advanced with a three-step low storage Runge-Kutta scheme. Optionally, for increased stability at low Reynolds numbers, at the price of higher computational demand, the diffusion term can be treated implicitly.

## Usage

### Input file

The input file `dns.in` sets the physical and computational parameters. In the `examples/` folder are examples of input files for several canonical flows. See `src/INFO_INPUT.md` for a detailed description of the input file.

For the *multi-block* implementation in branches [`multi_block`](https://github.com/p-costa/SNaC/tree/multi_block) and [`multi_block_fft`](https://github.com/p-costa/SNaC/tree/multi_block), the block files in `geo/block.???` setup the flow geometry.

Files `out1d.h90` and `out3d.h90` in `src/` set which data are written in 1- and 3-dimensional output files, respectively. *The code should be recompiled after editing out?d.h90 files*.

### Compilation

The code should be compiled in `src/`. The prerequisites are the following:

 * MPI
 * [HYPRE](https://github.com/hypre-space/hypre)
 * OpenMP (optional)
 * FFTW (optional, for branch [`multi_block_fft`](https://github.com/p-costa/SNaC/tree/multi_block))

The Makefile in `src/` should be modified in agreement to the installation paths of each library. Also, the following preprocessor options are available:

 * `-D_TIMING`           : wall-clock time per time step is computed
 * `-D_IMPDIFF`          : diffusion term of the N-S equations is integrated in time with an implicit discretization (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `-D_SINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double)
 * `-D_FFT_?`, with  `?` being `X`, `Y` or `Z`: will use FFTs to solve the Poisson equation in the direction in question. This option is only valid for the multi-block implementation in branch [`multi_block_fft`](https://github.com/p-costa/SNaC/tree/multi_block)).

Typing `make run` will compile the code and copy the executable `cans` and input file `dns.in` to a `run/` folder.

### Running the code

Run the executable with `mpirun` with a number of tasks and shared threads complying to what has been set in the input file `dns.in` (or in the `geo/block.???` files in case of multi-block). Data will be written by default in a folder named `data/`, which must be located where the executable is run.

### Visualizing field data

See [`src/INFO_VISU.md`](src/INFO_VISU.md).

## Notes

I appreciate any feedback that can improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.

Please read the `LICENSE` file.

## Contributors

Pedro Costa (p.simoes.costa@gmail.com)
