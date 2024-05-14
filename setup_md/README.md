## setup_md.py

This script and its associated files are used to set up GROMACS MD simulations.
In addition to setting up simple MD runs with standard molecules, the script
can load in non-standard molecules with custom forcefield parameters generated
by charmm-gui. The core script should be forcefield-agnostic; however, it has
only been tested for charmm36 and custom charmm-gui parameters written out in
GROMACS format. 

This is a command-line program run by invoking the script `setup_md.py`.
The following command would set up an MD simulation in the directory
`md-input` from the file `PROTEIN.pdb`.

```
setup_md.py md-input --standard PROTEIN.pdb
```

You can add molecules with custom parameters as well. The following would build
a simulation with `PROTEIN.pdb` and a non-standard molecule whose coordinates 
are in `pose1.pdb` and whose custom parameters are in `pose1_params`.
 
```
setup_md.py md-input --standard PROTEIN.pdb --custom pose1.pdb --params pose1_params/
```

There are a variety of ways to control the behavior of the script. To see
information about the arguments, type `setup_md.py --help`. 

### Installation
The script assumes it is in the `setup_md` directory, so you should keep the
directory intact rather than copying out just the `setup_md.py` script. I'd
recommend putting the `setup_md` directory in a convenient location and then
calling the script with its path. For example, if you placed the `setup_md`
directory in your home directory, you would call:

```
~/setup_md/setup_md.py
```

You can also add the directory to your PATH environment variable. For the 
example above, you could add the following line to your .bashrc file. 

```
export PATH=$HOME/setup_md/:$PATH
```

If you did this, you could then simply type `setup_md.py` to run the script.

### Customizing simulation files
The script will copy the entire contents of `setup_md/run_files` into the final
directory. You can thus put whatever combination of mdp files, batch files, and
scripts that you want to have in your final simulation. 

The script comes with srun files for the University of Oregon talapas cluster. 
These need to be updated to have the correct PIRG before use. (They might need 
more extensive modification for other clusters). 

### Forcefield
The script looks in the `setup_md` directory for the forcefield to use. As of 
this writing, this is `charmm36-jul2022.ff`. If you want to set a new default,
change the `DEFAULT_FF_DIR` variable at the top of the script and place your 
new forcefield into the `setup_md` directory. Alternatively, you can specify
the forcefield for a specific run with the `--ff-source-dir` flag.
