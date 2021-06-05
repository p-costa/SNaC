# how to visualize the output binary files form *SNaC*

in addition to the binary files for visualization, *SNaC* generates a log file that contains information about the saved data; this information is then used to generate a single `Xdmf` metadata file for visualizing field data as a time series.

the steps are as follows:

1. after the simulation has run, copy the contents of `utils/visualize_fields/write_xdmf_all.py` to the simulation `data` folder;
2. run the file with `python write_xdmf_all.py`.
3. load the generated Xdmf (`*.xmf`) file using paraview/visit or other visualization software.

## example: how to visualize the default binary output

### 3D fields

when running the script `write_xdmf_all.py` we get the following prompts:

~~~
 $ python write_xdmf_all.py
 Name of the log file written by SNaC (excluding the block-specific suffix) [log_visu_3d]:
 Name to be appended to the grid files to prevent overwriting []:
 Block # 001 log files parsed and grid files generated.
 Block # 002 log files parsed and grid files generated.
 Block # 003 log files parsed and grid files generated.
 Name of the output file [viewfld_DNS.xmf]:
~~~

* the first input is the name of the file that logged the saved data;
* the second is a name to append to the grid files that are generated, which should change for different log files to prevent conflicts;
* the last is the name of the visualization file.

by pressing <kbd>enter</kbd> after each prompt, the default values in the square brackets are assumed by the script; these correspond to the default steps required for visualizing 3D field data; the data can then be visualized with, e.g., `$ paraview viewfld_DNS.xmf`.

### checkpoint files

a similar script also located in `utils/visualize_fields/`, named `write_xdmf_restart.py`, can be used to generate metadata that allows to visualize the field data contained in saved checkpoint files.

### read binary data for post-processing

under `utils/read_bindary_data` are MATLAB and python scripts which may be used to load binary field data for post-processing; see `tests/lid_driven_cavity/test.py` for an example use of the python script.
