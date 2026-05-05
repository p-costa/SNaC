# SNaC Input Files

SNaC reads two Fortran namelist files from the run directory:

1. `dns.nml`, with global DNS parameters shared by all blocks.
2. `blocks.nml`, with geometry, decomposition, boundary conditions, and grid choices for each block.

The default files live in `src/`, and `make run` copies them to `run/`.

## `dns.nml`

Example:

```fortran
&dns
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
visci = 100.
nstep = 100, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
icheck = 10, iout0d = 10, iout1d = 200000, iout2d = 50000000, iout3d = 10, isave = 20
bforce(1:3) = 0., 0., 0.
gacc(1:3) = 0., 0., 0.
nscal = 0
/

&scalar
iniscal(:) = 'zer'
alphai(:) = 1.
beta = 0.
ssource(:) = 0.
is_sforced(:) = F
scalf(:) = 0.
/

&hypre
hypre_solver_i = 2, hypre_tol = 1.e-4, hypre_maxiter = 50
/
```

`cfl`, `dtmax`, and `dt_f` control the time step. If `dt_f < 0`, SNaC uses `min(cfl*dt_cfl,dtmax)`, where `dt_cfl` is the stability-limited time step computed by the code. If `dt_f > 0`, SNaC uses `dt_f` as a fixed time step.

`visci` is the inverse kinematic viscosity. SNaC computes the kinematic viscosity as `visci**(-1)`.

`nstep`, `time_max`, and `tw_max` define the maximum number of steps, physical time, and wall-clock time in hours. `stop_type(1:3)` enables those three criteria.

`restart` loads `data/fld_b_###_<field>.bin` checkpoint files. `is_overwrite_save` controls whether saves overwrite the current checkpoint or keep numbered checkpoint history. `nsaves_max` limits the number of numbered checkpoint files when overwrite is disabled.

`icheck`, `iout0d`, `iout1d`, `iout2d`, `iout3d`, and `isave` set the cadence for stability checks, history output, 1D profiles, 2D slices, 3D fields, and checkpoint saves.

`bforce(1:3)` is a constant body force density, often used as a pressure-gradient forcing.

`nscal` sets the number of scalar fields. When `nscal > 0`, the `&scalar` namelist is read. `alphai(:)` is the inverse scalar diffusivity, and `iniscal(:)` selects each scalar initial condition. Supported scalar initial conditions include `zer`, `uni`, `cou`, and `dhc`.

When SNaC is compiled with `BOUSSINESQ_BUOYANCY=1` in `build.conf` or on the `make` command line, the first scalar is active and contributes the Boussinesq acceleration `-gacc(:)*beta*s`. Without that CPP macro, no buoyancy source is compiled in.

The optional `&hypre` namelist controls the iterative Poisson solver. `hypre_solver_i` selects `1: SMG`, `2: PFMG`, `3: GMRES`, or `4: BiCGSTAB`.

## `blocks.nml`

Each block is indexed in the second or final namelist dimension. This example defines one block:

```fortran
&blocks
nblocks = 1
block_dims(1:3,1) = 1, 1, 1
block_ng(1:3,1) = 32, 32, 64
block_lmin(1:3,1) = 0., 0., 0.
block_lmax(1:3,1) = 0.5, 0.5, 1.
block_gt(1:3,1) = 0, 0, 0
block_gr(1:3,1) = 0., 0., 0.
block_cbcvel(0:1,1:3,1,1) = 'D', 'D', 'D', 'F', 'D', 'D'
block_cbcvel(0:1,1:3,2,1) = 'D', 'D', 'D', 'F', 'D', 'D'
block_cbcvel(0:1,1:3,3,1) = 'D', 'D', 'D', 'F', 'D', 'D'
block_cbcpre(0:1,1:3,1) = 'N', 'N', 'N', 'F', 'N', 'N'
block_bcvel(0:1,1:3,1,1) = 0., 0., 0., 2., 0., 0.
block_bcvel(0:1,1:3,2,1) = 0., 0., 1., 2., 0., 0.
block_bcvel(0:1,1:3,3,1) = 0., 0., 0., 2., 0., 0.
block_bcpre(0:1,1:3,1) = 0., 0., 0., 2., 0., 0.
block_cbcscal(0:1,1:3,1,1) = 'N', 'N', 'N', 'F', 'N', 'N'
block_bcscal(0:1,1:3,1,1) = 0., 0., 0., 2., 0., 0.
block_inflow_type(0:1,1:3,1) = 0, 0, 0, 0, 0, 0
block_inivel(1) = 'zer'
/
```

`block_dims(1:3,iblock)` is the MPI decomposition of block `iblock`.

`block_ng(1:3,iblock)` is the number of grid cells along each block direction.

`block_lmin` and `block_lmax` are the physical lower and upper corners of the block.

`block_gt` selects the SNaC mapping function used by `initgrid.f90`; `block_gr` is the corresponding stretching parameter. `gr = 0` gives a uniform grid.

The boundary-condition arrays use face order `x-`, `x+`, `y-`, `y+`, `z-`, `z+`. For velocity arrays, the third index is the velocity component: `1` for `u`, `2` for `v`, `3` for `w`.

Boundary-condition types are:

* `F`: friend block connection.
* `D`: Dirichlet.
* `N`: Neumann.

For an `F` boundary condition, the corresponding value in `block_bcvel` or `block_bcpre` is the connected block index. Periodic directions also use `F`; in a single-block periodic direction, both periodic faces point back to that same block.

If scalars are enabled, `block_cbcscal(0:1,1:3,is,iblock)` and `block_bcscal(0:1,1:3,is,iblock)` set the boundary conditions for scalar `is` on block `iblock`.

`block_inflow_type` enables the built-in inflow profile generator on Dirichlet velocity boundaries. A non-positive value enforces the constant boundary value.

`block_inivel` chooses the initial velocity field. Supported values include `zer`, `uni`, `cou`, `poi`, `log`, `hcp`, `hcl`, and `tgv`.

## External Axis Grids

After building the regular `block_gt/block_gr` grid, SNaC checks for optional binary axis grids in `grid/`. For block 1, the recognized files are:

```text
grid/grid_x_b_001.bin
grid/grid_y_b_001.bin
grid/grid_z_b_001.bin
```

When one of these files exists, SNaC loads that axis from disk. The binary layout matches `save_grid`: face coordinates, center coordinates, face spacings, and center spacings for the block's `ng` entries in that direction.

The GUI in `utils/snac_grid_generator` writes `blocks.nml` and these optional grid files. This lets a block mix SNaC-native mapping functions on some axes with generated multi-grading on others.

## Connectivity

Connected blocks must be congruent on shared faces. Neighboring subdomains must share the same boundary, and transverse `block_ng`, `block_dims`, and grid spacings must match on friend faces.

## Restart Files

New checkpoints are written as one raw binary file per field and block: `fld_b_###_u.bin`, `fld_b_###_v.bin`, `fld_b_###_w.bin`, `fld_b_###_p.bin`, and, for scalars, `fld_b_###_s_001.bin`, `fld_b_###_s_002.bin`, etc.
