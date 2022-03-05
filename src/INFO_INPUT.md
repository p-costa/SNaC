# about `dns.in` and `geo/block.{???}` files

Setting up a multi-block simulation involves two input files:

1. `dns.in`,    which sets up the governing parameters that are common to all blocks, *global parameters*.
2. `block.???`, files located in the `geo/` folder, which set up *block-specific parameters*.

---
---

Let us start with block-independent parameters in `dns.in`. Consider the following input file `dns.in` as example:

~~~
 .95 1.0e5                   ! cfl, dtmin
 1. 1. 100.                  ! uref, lref, rey
 100    100. 0.1             ! nstep,time_max,tw_max
 T F F                       ! stop_type(1:3)
 F T                         ! restart,is_overwrite_save
10 10 20 5000 10000 2000     ! icheck, iout0d, iout1d, iout2d, iout3d, isave
 0. 0. 0.                    ! bforce(1:3)
 4                           ! nthreadsmax
~~~

---
---

~~~
.95 1.0e5                ! cfl, dtmin
~~~

This line controls the simulation time step.

The time step is set to be equal to `min(cfl*dtmax,dtmin)`, i.e. the minimum value between `dtmin` and `cfl` times the maximum allowable time step `dtmax` (computed every `ickeck` time steps; see below).
`dtmin` is therefore used when a constant time step, smaller than `cfl*dtmax`, is required. If not, it should be set to a high value so that the time step is dynamically adjusted to `cfl*dtmax`.

---

~~~
1. 1. 100.              ! uref, lref, rey
~~~

This line defines the flow governing parameters.

`uref`, `lref` and `rey` are a reference **velocity scale**, **length scale**, and **Reynolds number** defining the problem. The fluid kinematic viscosity is computed form these quantities.

~~~
100000 100. 0.1          ! nstep, time_max, tw_max
T F F                    ! stop_type(1:3)
F T                      ! restart,is_overwrite_save
~~~

These lines set the simulation termination criteria and whether the simulation should be restarted from a checkpoint file.

`nstep` is the **total number of time steps**.

`time_max` is the **maximum physical time**.

`tw_max` is the **maximum total simulation wall-clock time**.

`stop_type` sets which criteria for terminating the simulation are to be used (more than one can be selected, and at least one of them must be `T`)

* `stop_type(1)`, if true (`T`), the simulation will terminate after `nstep` time steps have been simulated;
* `stop_type(2)`, if true (`T`), the simulation will terminate after `time_max` physical time units have been reached;
* `stop_type(3)`, if true (`T`), the simulation will terminate after `tw_max` simulation wall-clock time (in hours) has been reached;

a checkpoint file `fld.bin` will be saved before the simulation is terminated.

`restart`, if true, **restarts the simulation** from a previously saved checkpoint file, named `fld.bin`.

`is_overwrite_save`, if true, overwrites the checkpoint file `fld.bin` at every save; if false, a symbolic link is created which makes `fld.bin` point to the last checkpoint file with name `fld_???????.bin`. In the latter case, to restart a run from a different checkpoint one just has to point the file `fld.bin` to the right file, e.g.: ` ln -sf fld_0000100.bin fld.bin`.

`nsaves_max` limits the number of saved checkpoint files, if `is_over_write_save` is false; a value of `0` or any negative integer corresponds to no limit, and the code uses the file format described above; otherwise, only `nsaves_max` checkpoint files are saved, with the oldest save being overwritten when the number of saved checkpoints exceeds this threshold; in this case, files with a format `fld_????.bin` are saved (with `????` denoting the saved file number), with `fld.bin` pointing to the last checkpoint file as described above; moreover, a file `log_saves.out` records information about the time step number and physical time corresponding to each saved file number.

---

~~~
10 10 20 5000 10000 2000 ! icheck, iout0d, iout1d, iout2d, iout3d, isave
~~~

These lines set the frequency of time step checking and output:

* every `icheck` time steps **the new time step size** is computed according to the new stability criterion and cfl (above);
* every `iout0d` time steps **history files with global scalar variables** are appended; currently the forcing pressure gradient and time step history are reported;
* every `iout1d` time steps **1d profiles** are written (e.g. velocity and its moments) to a file;
* every `iout2d` time steps **2d slices of a 3d scalar field** are written to a file;
* every `iout3d` time steps **3d scalar fields** are written to a file;
* every `isave`  time steps a **checkpoint file** is written (`fld_???????.bin`), and a symbolic link for the restart file, `fld.bin`, will point to this last save so that, by default, the last saved checkpoint file is used to restart the simulation.

1d, 2d and 3d outputs can be tweaked modifying files `out?d.h90`, and re-compiling the source. See also `output.f90` for more details.

---

~~~
0. 0. 0.                 ! bforce(1:3)
~~~

`bforce`, is a constant **body force density term** in the direction in question (e.g. the negative of a constant pressure gradient) that can be added to the right-hand side of the momentum equation. The three values correspond to three domain directions. NOTE: in a pressure-driven wall-bounded flow, only one type of flow forcing should be selected (bulk velocity or pressure gradient). If the streamwise bulk velocity is forced (by setting the is_forced parameter below `T`), bforce should be zero, and vice-versa.

---

~~~
4                        ! nthreadsmax
~~~

These lines set the grid of computational subdomains and maximum number of threads.

`nthreadsmax ` is the **maximum number OpenMP threads**.

---
---

The geometry, boundary and initial conditions, domain decompositions is set in a series of block files located in the `geo/` folder. Listed below for the case of a L channel:

`geo/block.001`:
~~~
1 1 1                    ! dims(1:3)
1  1  1                  ! lo(1:3)
32 32 64                 ! hi(1:3)
0. 0. 0.                 ! lmin(1:3)
.5 .5 1.                 ! lmax(1:3)
0 0 0                    ! gt(1:3)
0. 0. 0.                 ! gr(1:3)
D D  D F  D D            ! cbcvel(0:1,1:3,1) [u BC type]
D D  D F  D D            ! cbcvel(0:1,1:3,2) [v BC type]
D D  D F  D D            ! cbcvel(0:1,1:3,3) [w BC type]
N N  N F  N N            ! cbcpre(0:1,1:3  ) [p BC type]
0. 0.  0. 2.  0. 0.      !  bcvel(0:1,1:3,1) [u BC value]
0. 0.  1. 2.  0. 0.      !  bcvel(0:1,1:3,2) [v BC value]
0. 0.  0. 2.  0. 0.      !  bcvel(0:1,1:3,3) [w BC value]
0. 0.  0. 2.  0. 0.      !  bcpre(0:1,1:3  ) [p BC value]
0 0  0 0  0 0            !  is_inflow(0:1,1:3)
zer                      ! inivel
1                        ! id
~~~
`geo/block.002`:
~~~
1 1 1                    ! dims(1:3)
1  33 1                  ! lo(1:3)
32 64 64                 ! hi(1:3)
0. .5 0.                 ! lmin(1:3)
.5 1. 1.                 ! lmax(1:3)
0 0 0                    ! gt(1:3)
0. 0. 0.                 ! gr(1:3)
D F  F D  D D            ! cbcvel(0:1,1:3,1) [u BC type]
D F  F D  D D            ! cbcvel(0:1,1:3,2) [v BC type]
D F  F D  D D            ! cbcvel(0:1,1:3,3) [w BC type]
N F  F N  N N            ! cbcpre(0:1,1:3  ) [p BC type]
0. 3.  1. 0.  0. 0.      !  bcvel(0:1,1:3,1) [u BC value]
0. 3.  1. 0.  0. 0.      !  bcvel(0:1,1:3,2) [v BC value]
0. 3.  1. 0.  0. 0.      !  bcvel(0:1,1:3,3) [w BC value]
0. 3.  1. 0.  0. 0.      !  bcpre(0:1,1:3  ) [p BC value]
0 0  0 0  0 0            !  is_inflow(0:1,1:3)
zer                      ! inivel
2                        ! id
~~~
`geo/block.003`:
~~~
2 1 1                    ! dims(1:3)
33  33 1                 ! lo(1:3)
128 64 64                ! hi(1:3)
.5 .5 0.                 ! lmin(1:3)
2.0 1. 1.                ! lmax(1:3)
0 0 0                    ! gt(1:3)
0. 0. 0.                 ! gr(1:3)
F N  D D  D D            ! cbcvel(0:1,1:3,1) [u BC type]
F N  D D  D D            ! cbcvel(0:1,1:3,2) [v BC type]
F N  D D  D D            ! cbcvel(0:1,1:3,3) [w BC type]
F D  N N  N N            ! cbcpre(0:1,1:3  ) [p BC type]
2. 0.  0. 0.  0. 0.      !  bcvel(0:1,1:3,1) [u BC value]
2. 0.  0. 0.  0. 0.      !  bcvel(0:1,1:3,2) [v BC value]
2. 0.  0. 0.  0. 0.      !  bcvel(0:1,1:3,3) [w BC value]
2. 0.  0. 0.  0. 0.      !  bcpre(0:1,1:3  ) [p BC value]
0 0  0 0  0 0            !  is_inflow(0:1,1:3)
zer                      ! inivel
3                        ! id
~~~

---
---

~~~
1 1 1                    ! dims(1:3)
1  1  1                  ! lo(1:3)
32 32 64                 ! hi(1:3)
0. 0. 0.                 ! lmin(1:3)
.5 .5 1.                 ! lmax(1:3)
0 0 0                    ! gt(1:3)
0. 0. 0.                 ! gr(1:3)
~~~
These lines set the domain decomposition and computational grid for each block.

`dims` is the **number of computational** subdomains in each direction.


`lo(1:3)` and `hi(1:3)` are the **coordinates of the lower and upper corners** of the block in question, **in index space**.

`lmin(1:3)` and `lmax(1:3)` are the **physical coordinates of the lower and upper corners** of the block in question.

`gt` sets the **type of grid stretching** in each direction for this block; it can take three values:

* `0` grid cells clustered in both ends of the domain (default for any value different than `1` and `2`)
* `1` grid cells clustered in the middle of the domain
* `2` grid cells clustered in the upper end of the domain

`gr` is the **grid stretching parameter** that tweaks the non-uniform grid in the third direction; zero `gr` implies a uniform grid. See `initgrid.f90` for more details.

---

~~~
D D  D F  D D            ! cbcvel(0:1,1:3,1) [u BC type]
D D  D F  D D            ! cbcvel(0:1,1:3,2) [v BC type]
D D  D F  D D            ! cbcvel(0:1,1:3,3) [w BC type]
N N  N F  N N            ! cbcpre(0:1,1:3  ) [p BC type]
0. 0.  0. 2.  0. 0.      !  bcvel(0:1,1:3,1) [u BC value]
0. 0.  1. 2.  0. 0.      !  bcvel(0:1,1:3,2) [v BC value]
0. 0.  0. 2.  0. 0.      !  bcvel(0:1,1:3,3) [w BC value]
0. 0.  0. 2.  0. 0.      !  bcpre(0:1,1:3  ) [p BC value]
~~~

These lines set the boundary conditions (BC) for the block in question.

The **type** (BC) for each field variable is set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1` where,

* `X0` `X1` set the type of BC the field variable for the **lower** and **upper** boundaries in `x`;
* `Y0` `Y1` set the type of BC the field variable for the **lower** and **upper** boundaries in `y`;
* `Z0` `Z1` set the type of BC the field variable for the **lower** and **upper** boundaries in `z`.

The four rows correspond to the three velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`.

The following options are available:

* `F` connectivity BC (*friend* of another block);
* `D` Dirichlet;
* `N` Neumann.

The **last four rows** follow the same logic, but now for the BC **values**. For a `F` BC, the BC value corresponds to the *friend* block this block is connected to, in the direction in question. A periodicity boundary condition is naturally set if the blocks are cyclically connected with `F` boundary conditions.

---

~~~
0 0  0 0  0 0            !  is_inflow(0:1,1:3)
~~~
These lines set an inflow boundary condition for the block in question. Right now, if `is_inflow > 0` *and* a Dirichlet BC for the velocity is employed in that direction, a Poiseuille-type inflow of the kind `vel(x1,x2) = velref*(vel1(x1)**p1)*(vel2(x2)**p2)`, where `vel1` and `vel2` are Poiseuille profiles with unit mean, `velref` the Dirichlet BC set for the velocity above, and the exponent `p1` (`p2`) is set to `0` if the direction of coordinate `x1` (`x2`) is periodic, and `1` otherwise. Finally, `is_inflow <= 0` and a Dirichlet BC for the velocity is employed in that direction, SNaC will simply prescribe that BC, i.e. enforce a constant inflow velocity.

---

~~~
zer                      ! inivel
~~~

This line sets the initial velocity field for the block in question.

`initvel` **chooses the initial velocity field**. The following options are available:

* `zer`: zero velocity field
* `uni`: uniform velocity field equal to `uref`                                     ; streamwise direction in `x`
* `cou`: plane Couette flow profile with symmetric wall velocities equal to `uref/2`; streamwise direction in `x`
* `poi`: plane Poiseuille flow profile with mean velocity `uref`                    ; streamwise direction in `x`
* `log`: logarithmic profile with mean velocity `uref`                              ; streamwise direction in `x`
* `hcp`: half channel with plane Poiseuille profile and mean velocity `uref`        ; streamwise direction in `x`
* `hcl`: half channel with logarithmic profile and mean velocity `uref`             ; streamwise direction in `x`
* `tgv`: Taylor-Green vortex

See `initflow.f90` for more details.

---

~~~
1                        ! id
~~~

Finally, this line defines the **ID of the block**. However, the line is currently not used. Instead, the ID of the block is determined from the extension of the block file. For the example shown here, block \#1 is recognized as the block defined in file `geo/block.001`, and ditto for blocks \#2 and \#3.

## A note on the connectivity between blocks

It is required that connected blocks **are congruent at the boundary**. This requirement holds not only for the definition of the geometry, but also for the domain decomposition and computational grid: **neighboring subdomains must share the same boundary, and the computational setup correspond to a regular rectilinear structured grid**. Hence, the value of `dims`, and parameters `gr` and `gt` must be set with this in mind. If these conditions are not met, the code will return a runtime error.
