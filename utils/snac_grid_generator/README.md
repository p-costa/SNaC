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
- achieved-grid diagnostics with clickable block/face navigation;
- previewed grid-congruence repair with per-axis authority locks;
- existing-case import from `blocks.nml`, project JSON, and binary grids;
- exact-rank MPI decomposition balancing in automatic, 1D, 2D, or 3D mode;
- array, mirror, align, and face-snap tools for multi-block selections;
- common SNaC velocity/pressure boundary presets;
- autosave recovery and project-schema migration;
- headless check, update, repair, and decomposition workflows;
- structured-grid validation and transactional case writing.

Python 3.8 or newer is supported. Install the runtime dependency from the
repository root with:

```sh
python3 -m pip install -r utils/snac_grid_generator/requirements.txt
```

Browser and lint tests use the separate test manifest:

```sh
python3 -m pip install -r utils/snac_grid_generator/requirements-test.txt
```

The grid-generator Python sources and tests have an enforced Flake8 baseline:

```sh
python3 -m flake8 utils/snac_grid_generator tests/snac_grid_generator
```

Run the GUI from the repository root:

```sh
python3 -m utils.snac_grid_generator.server
```

The server prints the local URL. The GUI writes files through the Python backend
to the output directory shown in the toolbar, by default:

```text
generated/snac_grid_case/
```

The server accepts API writes only from its own browser origin and binds to a
loopback address by default. A non-loopback `--host` requires the explicit
`--allow-remote` option; the server does not provide user authentication or TLS,
so remote mode is intended only for a trusted network. Import and output paths
may reference any location accessible to the user running the process.

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

Check messages that identify blocks or faces are buttons. Selecting one checks
all affected blocks, switches to the relevant axis, centers the viewport target,
and outlines affected faces in red or amber.

The selected axis can be copied to every face-connected block with the same
extent in that axis using `Apply to aligned blocks`. For example, a Y grid can
be propagated across a row of blocks connected through X faces. It is not
copied through a face normal to the selected axis, because those blocks occupy
different intervals of that axis.

`Repair` checks all connected tangential grids together and previews the changes
needed to make them congruent. It also matches normal spacing when an interface
jump exceeds the warning threshold. By default the selected block is
authoritative; the lock button beside an axis makes that axis authoritative
regardless of the selection. Conflicting locked grids are reported and left
untouched. A monotone target keeps its profile when one endpoint is constrained;
otherwise repair uses smooth explicit coordinates with the same cell count and
block extent. External-grid writing is enabled when those coordinates are needed.

The MPI panel finds a decomposition whose per-block rank counts add exactly to
the requested total. `Auto` compares permitted dimensionalities, while `1D`,
`2D`, and `3D` require that number of active partition axes. The X/Y/Z switches
can exclude directions, notably an FFT-synthesis direction, which SNaC requires
to remain undistributed. The optimizer balances cell load while favoring compact
local domains and lower communication area. The minimum number of cells along
each split direction and an optional maximum local partition aspect are project
settings. A maximum aspect of zero disables that hard limit. The result shows each block's dimensions, cells-per-rank
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
is normally filled by `Update`. No-slip, free-slip, moving-wall, inlet, and
outlet presets can be applied to the checked blocks; inferred friend faces are
left untouched.

The viewport can color the selected block's boundary faces: velocity `D` is
red, velocity `N` is blue, friend boundaries are green, and faces with mixed
velocity-component codes are amber. MPI partition planes use the same
base-plus-remainder decomposition as `initmpi`, so they show the actual local
subdomain cuts for the selected block.
The grid-line button toggles a sparse 3D view of the selected block's actual
face coordinates, including imported and repaired explicit grids.
After `Check` or `Update`, connected faces are outlined by their normal-spacing
ratio: green below 1.5, amber above 1.5, and red above the warning ratio of 3.

The checkboxes in the block list define a multi-block selection and the
highlighted block remains the primary source. Selected blocks can be duplicated
as an array, mirrored about a coordinate plane, or aligned to the primary block.
One checked target can also be snapped to a face of the primary block. Mirroring
reverses the corresponding grading as well as the block geometry.

The select, move, and resize viewport tools edit block bounds directly. A move
or resize invalidates the last structural check and refreshes the grid preview.
Undo and redo retain up to 100 coalesced project edits, including decomposition
diagnostics; reset can be undone, while opening another project starts a fresh
history. Project edits are autosaved in browser storage. On the next launch the
GUI offers to restore or discard the saved session.

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

Each export is fully staged and validated before publication. Existing
generator-owned files are backed up during publication and restored if any
replacement fails. The hidden `.snac_grid_generator_manifest.json` records file
ownership so later exports can remove obsolete block grids without deleting
unrelated solver output. A malformed manifest, symbolic-link target, or
unmanaged file collision stops export and leaves the existing case untouched.

The binary grid files follow `save_grid` in `src/initgrid.f90`: for each axis,
the file contains the face coordinate, center coordinate, face spacing, and
center spacing arrays for the block's `ng` entries. Export uses little-endian
double precision. Import accepts little- or big-endian single or double
precision and validates all four arrays against one another; invalid files are
rejected instead of importing plausible face coordinates from a corrupted
payload.

You can also export from a saved project JSON:

```sh
python3 -m utils.snac_grid_generator.cli generated/snac_grid_case/snac_grid_project.json -o generated/snac_grid_case
```

To open an existing case in the GUI, put its directory in `Output` and use the
folder-import button. `snac_grid_project.json` is preferred when present;
otherwise the importer reconstructs the project from `blocks.nml`. Any
`grid/grid_[xyz]_b_###.bin` files, or `data/` equivalents, replace the axis
definition with their exact face coordinates. If one copy is invalid and the
other is valid, the valid copy is imported with a warning. If both are valid but
differ, `grid/` remains authoritative and the difference is reported. The same
conversion is available headlessly:

```sh
python3 -m utils.snac_grid_generator.cli import path/to/case -o imported-project.json
```

The original export command remains valid; `export` may also be written
explicitly as its first argument.

The project operations are also available without the GUI. Add `--json` for
machine-readable diagnostics, and use `-o` on mutating commands to write the
resulting project:

```sh
python3 -m utils.snac_grid_generator.cli check project.json --json
python3 -m utils.snac_grid_generator.cli update project.json -o updated.json
python3 -m utils.snac_grid_generator.cli repair project.json -o repaired.json
python3 -m utils.snac_grid_generator.cli decompose project.json -o decomposed.json
python3 -m utils.snac_grid_generator.cli migrate old-project.json -o project.json
```

Geometry, grading, repair, validation, and decomposition are format-neutral
modules. `topology.py` owns face connectivity and deterministic global-index
reconstruction, independent of block IDs. SNaC's native `gt/gr` mappings and
binary layout live in `snac_grid.py`; SNaC boundary presets live in
`snac_bc.py`; and case import/export remains in the SNaC adapters. The browser
entry point coordinates focused modules for project history, API transport,
block geometry, scene helpers, and SNaC-specific constants. The GUI's pinned
Three.js r164 and Lucide 0.468.0 assets are vendored under `static/vendor`, so
the editor does not require a CDN connection. Update Three.js core and both
control modules together, retain each upstream license, record the new versions
here, and run the complete browser smoke test after an asset update.

## SNaC Reader Hook

SNaC still generates grids from `gt/gr` by default. If a generated file such as
`grid/grid_x_b_001.bin` exists, the updated solver loads that axis from disk after
the normal `initgrid` call. This allows mixed usage: a block can use native SNaC
mapping on one axis and externally generated multi-grading on another axis.

The GUI can infer friend (`F`) boundaries for full-face-congruent neighboring
blocks. Connected blocks must still obey SNaC's existing constraints: matching
face geometry, matching transverse grid lines, and matching transverse
`ng`/`dims`.
