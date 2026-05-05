# SNaC Multi-Block Grid Generator

This utility builds rectilinear multi-block SNaC cases from a small project
model. It includes:

- a local 3D browser GUI served by Python;
- SNaC `blocks.nml` writers;
- SNaC-compatible binary grid writers, one file per block and direction;
- block scalar boundary rows for any project-level `nscal`;
- SNaC `initgrid` mapping functions `0` through `7`;
- OpenFOAM-style multi-grading segments `(length weight, cell weight, expansion ratio)`;
- geometric, tanh, and erf profiles for ratio-based generated grids;
- one-sided single-ratio controls with lower/upper-end bias;
- min/max spacing controls with lower-end, upper-end, both-ends, and middle clustering;
- spacing calculators for the common blockMesh parameters;
- structured-grid validation before writing case files.

Run the GUI from the repository root:

```sh
python3 -m utils.snac_grid_generator.server
```

The server prints the local URL. The GUI writes files through the Python backend
to the output directory shown in the toolbar, by default:

```text
generated/snac_grid_case/
```

The reset button clears the project to zero blocks. Exporting that empty project
is allowed and writes `blocks.nml` plus `snac_grid_project.json`.

The toolbar `Check` action validates the block layout, friend boundary
connectivity, transverse grid compatibility, and MPI partition compatibility.
The `Update` action infers full-face friend (`F`) boundaries and propagates MPI
partition counts across connected blocks where SNaC requires the neighboring
subdomains to line up. Export runs the same structured-grid checks and refuses
to write files when they fail. Large normal-spacing jumps across block
interfaces are reported as warnings.

Periodicity is controlled in the geometry panel. A periodic direction is only
valid when every block has the same min/max extent in that direction; `Update`
then marks the two periodic faces of each block as self-connected `F`
boundaries.

Boundary rows expose the velocity, pressure, and scalar condition codes plus
their numeric values. For `F` faces, the value is the neighboring block id and
is normally filled by `Update`.

For an exported case, the important files are:

```text
blocks.nml
grid/grid_x_b_001.bin
grid/grid_y_b_001.bin
grid/grid_z_b_001.bin
snac_grid_project.json
data/geometry_b_001.out
```

The `Grid files` selector controls where external grid files are copied.
`grid/` writes only the startup files that SNaC looks for before falling back to
`gt/gr`. `grid/ + data/` also duplicates those grid files under `data/`, matching
the grid snapshots SNaC writes with saved output.

The binary grid files follow `save_grid` in `src/initgrid.f90`: for each axis,
the file contains the face coordinate, center coordinate, face spacing, and
center spacing arrays for the block's `ng` entries, in double precision.

You can also export from a saved project JSON:

```sh
python3 -m utils.snac_grid_generator.cli generated/snac_grid_case/snac_grid_project.json -o generated/snac_grid_case
```

## SNaC Reader Hook

SNaC still generates grids from `gt/gr` by default. If a generated file such as
`grid/grid_x_b_001.bin` exists, the updated solver loads that axis from disk after
the normal `initgrid` call. This allows mixed usage: a block can use native SNaC
mapping on one axis and externally generated multi-grading on another axis.

The GUI can infer friend (`F`) boundaries for full-face-congruent neighboring
blocks. Connected blocks must still obey SNaC's existing constraints: matching
face geometry, matching transverse grid lines, and matching transverse
`ng`/`dims`.
