# SNaC Multi-Block Grid Generator

This utility builds rectilinear multi-block SNaC cases from a small project
model. It includes:

- a local 3D browser GUI served by Python;
- SNaC `blocks.nml` writers;
- SNaC-compatible binary grid writers, one file per block and direction;
- block scalar boundary rows for any project-level `nscal`;
- SNaC `initgrid` mapping functions `0` through `7`;
- OpenFOAM-style multi-region grading with optional spacing continuity;
- geometric, tanh, and erf monotone and symmetric profiles;
- interchangeable cell-count, ratio, and endpoint-spacing controls;
- achieved-grid diagnostics and interface-spacing warnings;
- previewed grid-congruence repair with per-axis authority locks;
- exact-rank MPI decomposition balancing in automatic, 1D, 2D, or 3D mode;
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

The reset button clears the project to zero blocks. That state can be saved as
project JSON, but `Check` and `Write` reject it because SNaC requires at least
one block.

The toolbar `Check` action validates the block layout, friend boundary
connectivity, transverse grid compatibility, and MPI partition compatibility.
The `Update` action infers full-face friend (`F`) boundaries and propagates MPI
partition counts across connected blocks where SNaC requires the neighboring
subdomains to line up. Export runs the same structured-grid checks and refuses
to write files when they fail. Large normal-spacing jumps across block
interfaces are reported as warnings.

The selected axis can be copied to every face-connected block with the same
extent in that axis using `Apply to aligned blocks`. For example, a Y grid can
be propagated across a row of blocks connected through X faces. It is not
copied through a face normal to the selected axis, because those blocks occupy
different intervals of that axis.

`Repair` checks all connected tangential grids together and previews the changes
needed to make them congruent. By default the selected block is authoritative;
the lock button beside an axis makes that axis authoritative regardless of the
selection. Conflicting locked grids are reported and left untouched.

The MPI panel finds a decomposition whose per-block rank counts add exactly to
the requested total. `Auto` compares permitted dimensionalities, while `1D`,
`2D`, and `3D` require that number of active partition axes. The X/Y/Z switches
can exclude directions, notably an FFT-synthesis direction, which SNaC requires
to remain undistributed. The optimizer balances cell load while favoring compact
local domains and lower communication area. It keeps at least four cells along
each split direction. The result shows each block's dimensions, cells-per-rank
imbalance, minimum local edge, and worst local aspect ratio. If no exact solution
exists, the current dimensions remain unchanged and nearby feasible totals are
offered as buttons that can be balanced immediately. Setting the target to zero
disables this additional export check.

## Grading Modes

`Native SNaC` exposes the mapping functions implemented by `initgrid`. These
grids can be represented directly by `gt` and `gr` in `blocks.nml`.

`Monotone grading` creates a one-way geometric, tanh, or erf profile. The block
length is fixed by its geometry, and any two independent controls can define
the spacing:

- cell count;
- upper/lower endpoint-spacing ratio;
- geometric cell-to-cell ratio;
- lower endpoint spacing;
- upper endpoint spacing.

The cell-to-cell ratio is available only for the geometric profile. A ratio
below one clusters cells toward the upper end; a ratio above one clusters them
toward the lower end. When the chosen pair implies a non-integer cell count,
the closest grid is generated and both the ideal count and fit error are shown.

`Symmetric grading` accepts absolute minimum and maximum spacing targets and
clusters the minimum spacing at either both ends or the middle. Its cell count
is likewise fitted to an integer. Older projects using one-sided min/max modes
are migrated to equivalent monotone lower/upper spacing controls when opened.

`Multi-region grading` divides an axis into normalized length and cell-count
weights. Each segment has an end/start spacing ratio. Enabling the link on a
segment makes its first cell match the previous segment's last cell; Python
solves the segment ratio and reports the achieved start, end, cell allocation,
and interface jump.

The preview always reports achieved cell count, lower and upper spacing,
minimum and maximum spacing, upper/lower expansion, maximum adjacent-cell
growth, ideal cell count, and fit error. These achieved values are authoritative
for the binary grid that will be written.

Periodicity is controlled in the geometry panel. `Update` matches congruent
block faces on the global minimum and maximum boundaries and marks each pair as
friend (`F`) boundaries. Blocks between those outer faces retain their normal
internal connections, while physical boundaries around holes remain
unconnected. A single block spanning a periodic direction is self-connected.

Boundary rows expose the velocity, pressure, and scalar condition codes plus
their numeric values. For `F` faces, the value is the neighboring block id and
is normally filled by `Update`.

The viewport can color the selected block's boundary faces: velocity `D` is
red, velocity `N` is blue, friend boundaries are green, and faces with mixed
velocity-component codes are amber. MPI partition planes use the same
base-plus-remainder decomposition as `initmpi`, so they show the actual local
subdomain cuts for the selected block.

The select, move, and resize viewport tools edit block bounds directly. A move
or resize invalidates the last structural check and refreshes the grid preview.
Undo and redo retain up to 100 coalesced project edits, including decomposition
diagnostics; reset can be undone, while opening another project starts a fresh
history.

For an exported case, the important files are:

```text
blocks.nml
grid/grid_x_b_001.bin
grid/grid_y_b_001.bin
grid/grid_z_b_001.bin
snac_grid_project.json
data/geometry_b_001.out
```

`Write external grid` controls whether binary and text grid files are emitted.
The `Grid files` selector controls where those external grid files are copied.
`grid/` writes only the startup files that SNaC looks for before falling back to
`gt/gr`. `grid/ + data/` also duplicates those grid files under `data/`, matching
the grid snapshots SNaC writes with saved output.

Each export is staged before it replaces the case files. The hidden
`.snac_grid_generator_manifest.json` records files owned by the generator so a
later export can remove obsolete block grids without deleting unrelated solver
output.

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
