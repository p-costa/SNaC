## Synopsis

**SNaC** is [**CaNS**](https://github.com/p-costa/CaNS) spelled backwards, and is a multi-block code for massively parallel direct numerical simulations (DNS) of fluid flows. *SNaC* aims at combining the versatility of a multi-block DNS solver, with the FFT-based acceleration used in CaNS.

The solver is able to simulate the flow in any three-dimensional multi-block structured Cartesian grid. However, if the geometry has one homogeneous, 'extruded' direction with constant grid spacing, SNaC can use a very fast solver that exploits FFTs in this direction. This is SNaC's *warp drive* :rocket:, as it yields a huge speedup in wall-clock time per time step.

**Reference**

P. Costa. *A FFT-accelerated multi-block finite-difference solver for massively parallel simulations of incompressible flows*.
Comput. Phys. Commun. 271 : 108194 (2022) [[DOI:10.1016/j.cpc.2021.108194]](https://doi.org/10.1016/j.cpc.2021.108194) [[arXiv:2106.03583]](https://arxiv.org/pdf/2106.03583.pdf).

## News
[27/04/2026] The input files now use Fortran namelists, `dns.nml` and `blocks.nml`, and the optional grid-generator utility can write multi-block `blocks.nml` files plus binary axis grids. See [`docs/INFO_INPUT.md`](docs/INFO_INPUT.md) for more details.

## Features

Some features are:

 * Multi-block, three-dimensional parallelization
 * Hybrid MPI/OpenMP parallelization
 * FFT-based synthesis of the Poissonn equation along one direction
 * HYPRE library used to solve Poisson/Helmholtz equations
 * Passive-scalar transport, with optional Boussinesq buoyancy from the first scalar
 * Parallel I/O using MPI I/O
 * A different canonical flow can be simulated just by changing the input files

## Motivation

*SNaC* is meant to serve as a multi-block DNS code for fast, massively-parallel simulations of single-phase flows, and as a solid base solver on top of which more complex phenomena can be implemented, such as numerical methods for multiphase flows.

## Method

The fluid flow is solved with a standard second-order finite-difference/-volume pressure correction scheme. Time is advanced with a three-step low storage Runge-Kutta scheme. Optionally, for increased stability at low Reynolds numbers, at the price of higher computational demand, the diffusion term can be treated implicitly.

## Usage

### Input files

The input file `dns.nml` sets the physical and computational parameters, while `blocks.nml` sets block-specific geometry, grids, and boundary conditions. In the `examples/` folder are examples of input files for several canonical flows. See [`docs/INFO_INPUT.md`](docs/INFO_INPUT.md) for a detailed description of the input files.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2-, and 3-dimensional output files, respectively. *The code should be recompiled after editing out?d.h90 files*.

### Compilation

#### Prerequisites
The prerequisites for compiling SNaC are the following:

 * MPI
 * [*HYPRE*](https://github.com/hypre-space/hypre)
 * OpenMP (optional)
 * *FFTW* (optional, in case FFT acceleration is used)

#### In short
For most systems, SNaC can be compiled from the root directory with:

```bash
make
```

The root `Makefile` creates local configuration files from `configs/defaults/` at the first compilation. The generated `build.conf` file selects the compiler, optimization profile, and pre-processing options such as `IMPDIFF`, `BOUSSINESQ_BUOYANCY`, and `FFT_AXIS`; see [`docs/INFO_COMPILING.md`](docs/INFO_COMPILING.md) for details.

Typing `make run` will compile the code and copy the executable `snac`, `dns.nml`, `blocks.nml`, and any optional `grid/` folder to a `run/` folder.

### Running the code

Run the executable with `mpirun` with a number of tasks and shared threads complying with the decomposition set in `blocks.nml`. Data will be written by default in a folder named `data/`, which must be located where the executable is run.

### Grid generator

Run the optional multi-block grid editor from the repository root with:

```bash
python3 -m pip install .
snac-grid-generator
```

The launcher accepts the same `--host`, `--port`, `--no-open`, and
`--allow-remote` options as `python3 -m utils.snac_grid_generator.server`.

The editor can draw and transform single or selected groups of blocks, configure SNaC-native or generated axis grids, apply common boundary presets, infer friend and periodic boundaries, check structured-grid constraints, and write `blocks.nml` plus optional binary axis grids. Its local grid library saves reusable axis gradings and complete case templates, including rescaling explicit coordinates to a new block extent. It autosaves locally and runs entirely from vendored static assets.

`Repair` makes tangential grids congruent across connected faces and matches large normal-spacing jumps. The selected block supplies an unlocked grid; locking an axis makes that block authoritative, and conflicting locks are reported instead of overwritten. Proposed changes are shown before they are applied.

The MPI panel balances a prescribed total rank count exactly. `Auto`, `1D`, `2D`, and `3D` control how many partition axes may be used, while the X/Y/Z switches restrict the allowed directions. Disable an FFT-synthesis direction because SNaC requires one MPI partition along that axis. The result reports per-block dimensions, cells per rank, and imbalance; an impossible request leaves the current decomposition unchanged and suggests nearby feasible rank counts. A target of zero disables the rank-count constraint.

The balancing search also favors compact subdomains with lower communication area. The minimum cells per split and optional maximum local aspect are configurable hard limits. Its summary reports the minimum local edge and maximum local aspect ratio, and nearby feasible totals can be applied directly.

The editor can import an existing case from `blocks.nml` and recover exact external grids from `grid/` or `data/` binaries. Sparse 3D grid lines show actual face coordinates, while checked interfaces are colored by their spacing jump. The quality panel reports grid, interface, and rank-balance metrics and downloads deterministic JSON or Markdown reports. Put the case directory in `Output` and use the folder-import button, or run `python3 -m utils.snac_grid_generator.cli import path/to/case -o imported-project.json`. The CLI also provides `check`, `update`, `repair`, `decompose`, `migrate`, and `report` commands with machine-readable diagnostics.

### Visualizing field data

See [`docs/INFO_VISU.md`](docs/INFO_VISU.md).

## Notes

I appreciate any feedback that can improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.

Please read the `LICENSE` file.

## Contributors

Pedro Costa (p.simoes.costa@gmail.com)
