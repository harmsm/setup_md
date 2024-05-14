#!/usr/bin/env python
"""
Generate GROMACS simulation inputs from pdb files and charmm-gui outputs.
"""
__author__ = "Michael J. Harms"
__date__ = "2024-04-18"
__version__ = "1.0.1"
        
DEFAULT_FF_DIR = "charmm36-jul2022.ff"
DEFAULT_RUN_FILES_DIR = "run_files"

import re
import glob
import os
import subprocess
import shutil
import sys
import argparse

class ITPManager:
    """
    Manage the content of ITP files. 
    """

    def _parse_simple(self,line):
        """
        Template parser. Single value key, values are all other columns.
        """
        
        col = line.split()
        key = col[0]
        values = tuple(col[1:])

        return key, values
    
    def _parse_pair_and_sort(self,line):
        """
        Template parser. Two value key, sorted. Values are all other columns.
        """

        col = line.split()
        
        key = line.split()[0:2]
        key.sort()
        key = tuple(key)
        
        values = tuple(col[2:])
        
        return key, values
    
    def _parse_on_number(self,line,number=1):
        """
        Template parser. Key is first col[:number] columns, value us col[number:].
        """

        col = line.split()
        key = tuple(col[:number])
        values = tuple(col[number:])

        return key, values

    def _parse_defaults(self,line):
        """
        Extract the key and values for an entry under a [ defaults ] field.
        """
        
        return self._parse_simple(line)

    
    def _parse_atomtypes(self,line):
        """
        Extract the key and values for an entry under an [ atomtypes ] field.
        """

        return self._parse_simple(line)
    
    def _parse_nonbond_params(self,line):
        """
        Extract the key and values for an entry under a [ nonbond_params ] field.
        """
        
        return self._parse_pair_and_sort(line)

    def _parse_bondtypes(self,line):
        """
        Extract the key and values for an entry under a [ bondtypes ] field.
        """
        
        return self._parse_pair_and_sort(line)

    def _parse_pairtypes(self,line):
        """
        Extract the key and values for an entry under a [ pairtypes ] field.
        """
        
        return self._parse_pair_and_sort(line)
    
    def _parse_angletypes(self,line):
        """
        Extract the key and values for an entry under an [ angletypes ] field.
        """
        
        return self._parse_on_number(line,number=3)
    
    def _parse_dihedraltypes(self,line):
        """
        Extract the key and values for an entry under a [ dihedraltypes ] field.
        """
        
        cols = line.split()

        key = cols[:4]
        values = cols[4:]
        if len(cols) == 8:
            key.append(cols[7])
            values = cols[4:7]

        return tuple(key), tuple(values)
    
    def _parse_cmaptypes(self,line):
        """
        Extract the key and values for an entry under a [ cmaptypes ] field.
        """

        return self._parse_on_number(line,number=5)

    def _parse_moleculetype(self,line):
        """
        Extract the key and values for an entry under a [ moleculetype ] field.
        """

        return self._parse_simple(line)

    def _parse_atoms(self,line):
        """
        Extract the key and values for an entry under an [ atoms ] field.
        """

        cols = line.split()
        key = (cols[3],cols[4])
        values = tuple(cols)

        return key, values
    
    def _parse_bonds(self,line):
        """
        Extract the key and values for an entry under a [ bonds ] field.
        """

        return self._parse_on_number(line,number=2)
    
    def _parse_pairs(self,line):
        """
        Extract the key and values for an entry under a [ pairs ] field.
        """

        return self._parse_on_number(line,number=2)
    
    def _parse_angles(self,line):
        """
        Extract the key and values for an entry under an [ angles ] field.
        """

        return self._parse_on_number(line,number=3)
    
    def _parse_dihedrals(self,line):
        """
        Extract the key and values for an entry under a [ dihedrals ] field.
        """

        return self._parse_on_number(line,number=4)

    def _parse_cmap(self,line):
        """
        Extract the key and values for an entry under a [ cmap ] field.
        """

        return self._parse_on_number(line,number=4)

    def _parse_dihedral_restraints(self,line):
        """
        Extract the key and values for an entry under a [ dihedral_restraints ]
        field.
        """

        return self._parse_on_number(line,number=4)
    
    def _parse_position_restraints(self,line):
        """
        Extract the key and values for an entry under a [ position_restraints ]
        field.
        """

        return self._parse_simple(line)
    
    def _parse_system(self,line):
        """
        Extract the key and values for an entry under a [ position_restraints ]
        field.
        """

        return self._parse_generic(line)
    
    def _parse_molecules(self,line):
        """
        Extract the key and values for an entry under a [ position_restraints ]
        field.
        """

        return self._parse_simple(line)
        
        
    def _parse_generic(self,line):
        """
        Extract the key and values for an entry under an unknown field. Return 
        the whole line as both the key and value. 
        """

        out = line.split()
        return tuple(out), tuple(out)

    def __init__(self,itp_file):

        # Store itp file
        self._itp_file = itp_file
        
        self._entry_types = ["defaults",
                             "atomtypes",
                             "nonbond_params",
                             "bondtypes",
                             "pairtypes",
                             "angletypes",
                             "dihedraltypes",
                             "cmaptypes",
                             "moleculetype",
                             "atoms",
                             "bonds",
                             "pairs",
                             "angles",
                             "dihedrals",
                             "cmap",
                             "dihedral_restraints",
                             "position_restraints",
                             "system",
                             "molecules"]

        # Build list of functions for parsing fields
        self._entry_parsers = {}
        for m in self._entry_types:
            self._entry_parsers[m] = getattr(self,f"_parse_{m}")

        # Fields that should be merged so they are only present in one itp file
        self._merge_fields = ["defaults",
                              "atomtypes",
                              "nonbond_params",
                              "bondtypes",
                              "pairtypes",
                              "angletypes",
                              "dihedraltypes",
                              "cmaptypes"]

        self._load_lines()
        self._load_fields()

    
    def _load_lines(self):
        """
        Load the lines of the ITP file. Creates _real_lines (holding lines with
        line breaks merged) and _raw_lines (dictionary with actual lines keyed
        by real line number). """
    
        with open(self._itp_file) as f:
        
            # raw_lines is a dictionary that maps real line numbers to lists of
            # raw_lines. For example, if we merge 5 lines that end with \ into
            # a single real line, the merged (real) line number will be a key
            # to a list holding the 5 raw lines in raw_lines. If we want to 
            # delete that real line (say, a cmaptypes entry spread across dozens
            # of raw lines), we simply ignore the real line key. 
            
            real_lines = []
            raw_lines = {len(real_lines):[]}
            
            line_queue = []
            for raw_line in f:
        
                raw_lines[len(real_lines)].append(raw_line)
                
                # This block of code merges lines that end in \. We take everything 
                # before ';' (comment) on each raw_line and then strip leading/trailing
                # spaces. If stripped_line ends with \, it is placed into the line_queue 
                # and we move to the next raw_line. If stripped_line does not end with
                # \, we append to line_queue, then merge the line_queue (with however 
                # many previous lines where in it). The merged line_queue is the 
                # real_line holding the complete input. It can be mapped back to the 
                stripped_line = raw_line.split(";")[0].strip()
                if len(stripped_line) > 0 and stripped_line[-1] == "\\":
                    line_queue.append(stripped_line[:-1])
                    continue
                else:
                    line_queue.append(stripped_line)
        
                real_line = " ".join(line_queue)
                line_queue = []
    
                real_lines.append(real_line)
                raw_lines[len(real_lines)] = []
                
        self._real_lines = real_lines
        self._raw_lines = raw_lines
        self._real_line_index = range(len(raw_lines))

    def _load_fields(self):
        """
        Load the data from an itp file into self._fields (keying field and
        entry to value) and self._line_numbers (keying field and entry to 
        real line number).
        """

        self._fields = {}
        self._line_numbers = {}

        # Lines with [ xxx ] (field_lines, nothing, or starting with #)
        self._field_lines = {} 
        self._blank_lines = []
        self._hash_lines = []
        
        field = None
        for line_counter, line in enumerate(self._real_lines):
            
            # Empty line, move on
            if len(line.strip()) == 0:
                self._blank_lines.append(line_counter)
                continue

            if line.startswith("#"):
                self._hash_lines.append(line_counter)
                continue
        
            # Look for something matching [ x ]. 'x' is the field
            if re.match("\\[ .*? \\]",line):
                
                field = line.strip()[1:-1].strip()
                if field not in self._fields:
                    self._fields[field] = {}
                    self._line_numbers[field] = {}
                    self._field_lines[field] = [line_counter]

                else:
                    self._field_lines[field].append(line_counter)
                    
                continue
                
            if field is None:
                err = "mangled ITP file? We should not have gotten here without\n"
                err += "figuring out the field.\n"
                raise ValueError(err)

            # Figure out how to parse this field
            if field in self._entry_parsers:
                entry, value = self._entry_parsers[field](line)
            else:
                entry, value = self._parse_generic(line)
            
            # Record value and line from which we got this entry. 
            if entry not in self._fields[field]:
                self._fields[field][entry] = [value]
                self._line_numbers[field][entry] = [line_counter]
            else:
                if issubclass(type(self._fields[field][entry]),tuple):
                    self._fields[field][entry] = [self._fields[field][entry]]
                    self._line_numbers[field][entry] = [self._line_numbers[field][entry]]

                self._fields[field][entry].append(value)
                self._line_numbers[field][entry].append(line_counter)
            
    @property
    def fields(self):
        """
        Dictionary where fields and entries key to lists of values. 
        """
        return self._fields
    
    @property
    def line_numbers(self):
        """
        Dictionary where field and entries key to lists of real line numbers.
        """
        return self._line_numbers

    def entries_in_other(self,other_im):
        """
        Get list of lines in self that that repeat field/entry pairs in some 
        other ITPManager instance.

        Parameters
        ----------
        other_im : ITPManager
            Initialized ITPManager instance against which to compare self.

        Returns
        -------
        shared_lines : list
            list of real line numbers in self that correspond to field/entry 
            pairs somewhere in the other ITP file. 
        """

        if not issubclass(type(other_im),type(self)):
            err = f"other_im should be an instance of {type(self)}"
            raise ValueError(err)

        shared_lines = []
        for field in self.fields:

            # Only merge fields if in self._merge_fields
            if field not in self._merge_fields:
                continue

            if field in other_im.fields:
                entries_here = set(self.fields[field])
                entries_there = set(other_im.fields[field])
                shared_entries = entries_here.intersection(entries_there)
                
                for entry in shared_entries:
                    shared_lines.extend(self.line_numbers[field][entry])
                        
        shared_lines.sort()
        
        return shared_lines

    def mask_existing_entries(self,other_itp_file):
        """
        Mask any lines that correspond to field/entry pairs found in some other
        itp file. Removes duplicated lines from self._real_line_index so, when
        we call self.write_itp, those entries are not written out. 
        
        Parameters
        ----------
        other_itp_file : str
            some other ITP file that could have duplicate entries relative to 
            the current one
        """

        other_im = ITPManager(other_itp_file)
        shared_lines = self.entries_in_other(other_im)
        
        new_line_numbers = set(self._real_line_index) - set(shared_lines)
        new_line_numbers = list(new_line_numbers)
        new_line_numbers.sort()

        self._real_line_index = new_line_numbers

    def mask_field(self,field):
        """
        Mask an entire field, including the [ field ] line. This removes lines
        from self._real_line_index so those lines are not written when we call
        self._write_itp.
        
        Parameters
        ----------
        field : str
            field to mask. must have been seen in the intial itp file. 
        """

        lines_to_remove = []
        for entry in self.line_numbers[field]:
            lines_to_remove.extend(self.line_numbers[field][entry])
        lines_to_remove = set(lines_to_remove)

        new_line_numbers = set(self._real_line_index) - lines_to_remove
        new_line_numbers = list(new_line_numbers)
        for line in self._field_lines[field]:
            new_line_numbers.remove(line)
        new_line_numbers.sort()
        self._real_line_index = new_line_numbers

    def write_itp(self,
                  itp_file,
                  masked=True):
        """
        Write out the ITP file. 
        
        Parameters
        ----------
        itp_file : str
            file to write to. (does not check for existing file and will
            overwrite)
        masked : bool, default=True
            if masked is True, we will not write lines previously masked via 
            mask_existing_entries and mask_field. If False, all lines will be
            written out, even if masked previously. 
        """

        if masked:
            line_numbers = self._real_line_index
        else:
            line_numbers = range(len(self._real_lines))

        out = []
        for i in line_numbers:
            out.extend(self._raw_lines[i])

        with open(itp_file,'w') as f:
            f.writelines(out)

def _insert_new_include(what_to_include,
                        itp_file,
                        insert_index=0,
                        insert_before=False):
    """
    Insert a new #include "what_to_include" directive into an itp file.

    Parameters
    ----------
    what_to_include : str
        file to include (via new line like '#include "what_to_include"')
    itp_file : str
        itp file to edit with the new #include
    insert_index : int, default=0
        insert the #include statement adjacent to the insert_index-th existing
        include. insert_index=0 means insert adjacent to the first #include. if
        insert_index > than the number of #include already there, place at end
    insert_before : bool, default=False
        If True, insert the #include before the #include at insert_index. The 
        default is to insert after.
    """

    # Whether to insert before or after #include found
    if insert_before:
        offset = 0
    else:
        offset = 1

    # Read in the topology file
    with open(itp_file,"r") as f:
        itp_lines = f.readlines()

    # Find lines that start with "#include". 
    include_lines = [i for i in range(len(itp_lines))
                     if itp_lines[i].startswith("#include")]
    
    # Figure out what index we are going to use for the insertion
    if len(include_lines) > insert_index and len(include_lines) > 0:
        idx = include_lines[insert_index] + offset
    else:
        idx = -1

    # Insert the new line
    itp_lines.insert(idx,f"\n#include \"{what_to_include}\"\n")

    # Write the altered file
    with open(itp_file,"w") as f:
        f.writelines(itp_lines)

def _get_included_files(some_itp,
                        seen_itp_files=None,
                        base_path=None):
    """
    Recursively search an input itp or top file for all files included using an
    #include directive.
    
    Parameters
    ----------
    some_itp : str
        itp or top file to search
    seen_itp_files : list, optional
        list of previously visited included files. This is used in the 
        recursion; it is usually not necessary for a user to set. 
    base_path : str, optional 
        directory name holding itp file that referenced the current file. This
        is used in the recursion; it is usually not necessary for a user to set.

    Returns
    -------
    seen_ipt_files : list
        list of unique itp files, including the initial file
    """

    # initial call. populate seen_itp_files and base_path for recursion. 
    if seen_itp_files is None:
        seen_itp_files = [some_itp]

    if base_path is None:
        base_path = os.path.relpath(os.path.dirname(some_itp))

    # Go through the itp file
    with open(some_itp) as f:
        
        for line in f:

            # Parse #include line
            if line.startswith("#include"):
                col = line.split()
                if len(col) == 2:
                    
                    # Get itp file (with full path given base_path)
                    itp_file = col[1].strip("\"")
                    dirname = os.path.dirname(itp_file)
                    basename = os.path.basename(itp_file)
                    
                    itp_file = os.path.join(base_path,dirname,basename)
                    itp_file =  os.path.relpath(itp_file)
                    
                    # If we haven't already seen itp_file, call function 
                    # recursively. 
                    if itp_file not in seen_itp_files:
                        seen_itp_files.append(itp_file)
                        seen_itp_files = _get_included_files(some_itp=itp_file,
                                                             seen_itp_files=seen_itp_files,
                                                             base_path=os.path.join(base_path,dirname))

    return seen_itp_files

def _look_for_pattern_in_files(text_search_pattern,
                               file_search_pattern="*.mdp",
                               search_dir="."):
    """
    Look for a text_search_pattern in files that match file_search_pattern
    in search_dir. (Recursive grep call.) 

    Parameters
    ----------
    text_search_pattern : str
        look in files for this pattern
    file_search_pattern : str
        look for file names that match this pattern
    search_dir : str, default='.'
        recursively look for files in this directory

    Returns
    -------
    hits : list
        list of tuples holding hits. tuples have form (filename_with_match,
                                                       line_number_with_match,
                                                       re.Match_object)
    """
    
    search_for = re.compile(text_search_pattern)

    search_in = os.path.join(search_dir,"**",file_search_pattern)
    matching_files = glob.glob(search_in,recursive=True)
    
    # Go through files and search for hits, building a list of tuples of
    # matches. Each tuple is (filename,line number of match,re_match_object). 
    hits = []
    for m in matching_files:
        with open(m) as f:
            for i, line in enumerate(f.readlines()):
                hit = search_for.search(line)
                if hit:
                    hits.append((m,i,hit))

    return hits

def _run_cmd(cmd_list,
             base=None,
             run_in_dir=None,
             input=None,
             print_stderr=False,
             print_stdout=False):
    """
    Run a bash command, doing error checking. 

    Parameters
    ----------    
    cmd_list : list
        bash command, with each argument an individual entry in the list.
    base : str, optional
        base to write out logs.  If none, do not write. 
    run_in_dir : str, optional
        directory in which to run the command. if none, run in current
        directory
    input : bytes, optional
        pipe input into bash call
    print_stderr : bool, default=False
        print and flush standard error
    print_stdout : bool, default=False
        print and flush standard output
    """

    # Go into specified working directory
    cwd = os.getcwd()
    if run_in_dir is not None:
        os.chdir(run_in_dir)
    
    # Run subprocess
    s = subprocess.run(cmd_list,capture_output=True,input=input)
    
    # Check for success
    if s.returncode != 0:
        err = ["Command failed.  Command was:\n\n"]
        err.append(" ".join(cmd_list))
        err.append("\n\n")
        err.append("standard error was:\n\n")

        stderr = s.stderr.decode("utf-8")
        
        err.append(stderr)
        
        os.chdir(cwd)
        
        raise RuntimeError("".join(err))
    
    # Write out logs
    if base:
        f = open(f"{base}.stdout","w")
        f.write(s.stdout.decode("utf-8"))
        f.close()
        
        f = open(f"{base}.stderr","w")
        f.write(s.stderr.decode("utf-8"))
        f.close()
   
    if print_stdout:
        print(s.stdout.decode("utf-8"),flush=True)
    
    if print_stderr:
        print(s.stderr.decode("utf-8"),flush=True)
 
    os.chdir(cwd)


def _offset_pdb_residues(pdb_file):
    """
    Offset the chains in a pdb file by adding 1000 to indicate chain. This is 
    helpful as chain labels can be lost during a simulation.
    
    Parameters
    ----------
    pdb_file : str
        pdb file to offset

    Returns
    -------
    out_file : str
        path to pdb file containing offset residues
    """

    # Get unique chains in the pdb file
    all_chains = []
    with open(pdb_file) as f:
        
        for line in f:
            if line[:6] in ["HETATM","ATOM  "]:
                chain = line[21]
                all_chains.append(chain)

    all_chains = [c for c in list(set(all_chains)) if c.strip() != ""]
    all_chains.sort()

    # No offset needed -- we only have one chain
    if len(all_chains) == 1:
        return pdb_file

    # Calculate offsets for chains
    offsets = {}
    current_offset = 0
    for c in all_chains:
        offsets[c] = current_offset
        current_offset += 1000

    # Do offsets
    out_lines = []
    with open(pdb_file) as f:
        for line in f:
            if line[:6] not in ["HETATM","ATOM  "]:
                out_lines.append(line)
                continue

            # Get current chain
            chain = line[21].strip()
            if chain == "":
                out_lines.append(line)
                continue

            # Offset residue. Check to see if residue overflows the residue field.
            resid = int(line[22:26])
            resid = resid + offsets[chain]
            if resid > 9999:
                err = "pdb resid count exceeded 9999. Try changing your residue\n"
                err += "numbers manually before running this script or run with\n"
                err += "no resid_offset applied.\n"
                raise ValueError(err)

            # Record updated line
            out_lines.append(f"{line[:22]}{resid:-4d}{line[26:]}")

    # Figure out new file name
    split_file = pdb_file.split(".")
    out_file = f"{'.'.join(split_file[:-1])}_offset.{split_file[-1]}"

    # Write to file
    with open(out_file,'w') as f:
        f.writelines(out_lines)

    return out_file


def initialize_directory(standard_pdb_file,
                         working_dir,
                         water_model,
                         ff_source_dir,
                         use_resid_offset,
                         prefix):
    """
    Create an initial topol.top and system.gro file. 

    Parameters
    ----------
    standard_pdb_file : str or None
        pdb file with only standard atoms. If defined, run this through pdb2gmx.
        If None, build an initial empty gro file and topology that includes the
        forcefield, water, and ions.
    working_dir : str
        directory to initialize
    water_model : str
        water model to use
    ff_source_dir : str
        forcefield directory name within working. we assume it ends in .ff
    use_resid_offset : bool
        whether or not to offset residue numbers by chain in standard_pdb_file 
    prefix : str
        pre-pend this to all log files so we can tell what run this is
        associated with
    """
    
    # Decide on the water model
    water_model_files = {"spc":"spc.itp",
                         "spce":"spce.itp",
                         "tip3p":"tip3p.itp",
                         "tip4p":"tip4p.itp",
                         "tip5p":"tip5p.itp"}

    # Throw error if water model is not recognized
    if water_model not in water_model_files:
        err = f"water_model {water_model} not recognized.\n"
        err += "should be one of:\n"
        for m in water_model_files:
            err += f"    {m}\n"
        raise ValueError(err)

    water_model_file = water_model_files[water_model]

    # .gro and .top files for the system. 
    system_gro = os.path.join(working_dir,"system.gro")
    system_topol = os.path.join(working_dir,"topol.top")

    # If there is a standard pdb file, run pdb2gmx
    if standard_pdb_file is not None:

        print(f"Loading {standard_pdb_file}",flush=True)
        
        tmp_standard_pdb = os.path.join(working_dir,f"{prefix}_standard.pdb")
        shutil.copy(standard_pdb_file,tmp_standard_pdb)

        if use_resid_offset:
            tmp_standard_pdb = _offset_pdb_residues(tmp_standard_pdb)
        
        local_ff = ff_source_dir.split(".ff")[0].strip()

        _run_cmd(["gmx","pdb2gmx",
                  "-f",os.path.split(tmp_standard_pdb)[-1],
                  "-o","system.gro",
                  "-p","topol.top",
                  "-ignh", 
                  "-ff",local_ff,
                  "-water",water_model,
                  "-cmap"],
                  base=f"{prefix}_pdb2gmx",
                  run_in_dir=working_dir,
                  print_stderr=True)
        
    # Otherwise, create a starting empty gro and topology file
    else:
        
        # Initial .gro with no atoms
        with open(system_gro,"w") as f:
            f.write(f"INITIAL GRO FILE\n{0:5d}\n{0:10.5f}{0:10.5f}{0:10.5f}")
        
        # Initial .top with core forcefield, water model, and ions. 
        with open(system_topol,"w") as f:
            f.write(f"#include \"./{ff_source_dir}/forcefield.itp\"\n")
            f.write(f"#include \"./{ff_source_dir}/{water_model_file}\"\n")
            f.write("\n\n")
            f.write("#ifdef POSRES_WATER\n")
            f.write("[ position_restraints ]\n")
            f.write("1    1       1000       1000       1000\n")
            f.write("#endif\n\n")
            f.write(f"#include \"./{ff_source_dir}/ions.itp\"\n")
            f.write("\n\n")
            f.write("[ system ]\n; Name\ngmx\n\n[ molecules ]\n\n")

def create_pre_solvate_ndx(working_dir,
                           prefix):

    cmd_file = f"{prefix}_ndx.command"
    with open(os.path.join(working_dir,cmd_file),"w") as f:
        f.write("name 0 user_input\n")
        f.write("keep 0\n")
        f.write("q\n")

    pipe_cmd = ["name 0 pre_solvate","keep 0","q",""]
    pipe_cmd = "\n".join(pipe_cmd)
    
    _run_cmd(["gmx","make_ndx",
              "-f","system.gro",
              "-o","pre_solvate.ndx"],
             base=f"{prefix}_make_ndx",
             run_in_dir=working_dir,
             input=pipe_cmd.encode('utf-8'),
             print_stderr=False)

def solvate_and_ionize(working_dir,
                       tmp_base,
                       box_size=1.0,
                       ion_conc=0.1,
                       pcharge=1,
                       ncharge=-1,
                       pname="NA",
                       nname="CL"):
    """
    Solvate and ionize a structure. This assumes that "working_dir" has 
    (minimally) system.gro and topol.top. Presumably the .top file has also 
    included all .itp files necessary to specify the system. 
    
    Parameters
    ----------
    working_dir : str
        directory in which to do the calculation
    tmp_base : str
        pre-pend this to all log files so we can tell what run this is
        associated with
    box_size : float, default=1
        box size in nm
    ion_conc : float, default=0.1
        ion concentration in molar
    pcharge : int, default=1
        charge of positive ions
    ncharge : int, default=-21
        charge of negative ions
    pname : str, default="NA"
        name of positive ions used to neutralize system
    nname : str, default="CL"
        name of negative ions used to neutralize system
    """

    print("Create box",flush=True)
    _run_cmd(["gmx","editconf",
             "-f","system.gro",
             "-o","system.gro",
             "-c",
             "-d","{:.1f}".format(box_size),
             "-bt","cubic"],
             base=f"{tmp_base}_editconf",
             run_in_dir=working_dir)
    
    print("Solvating",flush=True)
    _run_cmd(["gmx","solvate",
             "-cp","system.gro",
             "-cs","spc216.gro",
             "-o","system.gro",
             "-p","topol.top"],
             base=f"{tmp_base}_solvate",
             run_in_dir=working_dir)
    
    print("Ionizing",flush=True)
    _run_cmd(["gmx","grompp",
             "-f","ions.mdp",
             "-c","system.gro",
             "-p","topol.top",
             "-o","ions.tpr"],
             base=f"{tmp_base}_grompp_ion",
             run_in_dir=working_dir)
            
    _run_cmd(["gmx","genion",
             "-s","ions.tpr",
             "-o","system.gro",
             "-p","topol.top",
             "-pq",f"{int(round(pcharge,0)):d}",
             "-nq",f"{int(round(ncharge,0)):d}",
             "-pname",pname,
             "-nname",nname,
             "-neutral",
             "-conc","{:.3}".format(ion_conc)],
             base=f"{tmp_base}_genion",
             run_in_dir=working_dir,
             input=b"SOL")

def import_custom_params(gromacs_param_dir,
                         working_dir,
                         prefix):
    """
    Import a `gromacs` directory written out by charmm-gui into the working
    directory. Creates an output directory with four files: molecule.itp,
    forcefield.itp, molecule.gro, and posres.txt

    Parameters
    ----------
    gromacs_param_dir : str
        gromacs directory from a charmm-gui forcefield conversion run
    working_dir : str
        working directory for the calculation
    prefix : str
        pre-pend this to all log files so we can tell what run this is
        associated with

    Returns
    -------
    output_directory : str
        directory in which we placed this imported set of parameters
    """

    # -------------------------------------------------------------------------
    # Make sure the input directory has appropriate files and structure. 

    # Look for gro file for custom molecule
    gro_files = glob.glob(os.path.join(gromacs_param_dir,"*.gro"))
    if len(gro_files) != 1:
        err = f"There should be a single .gro file in {gromacs_param_dir}\n"
        raise ValueError(err)
    gro_file = gro_files[0]

    # Look for itp file for custom molecule
    itp_files = glob.glob(os.path.join(gromacs_param_dir,"toppar","*.itp"))
    itp_files = [itp for itp in itp_files if os.path.split(itp)[-1] != "forcefield.itp"]
    if len(itp_files) != 1:
        err = f"There should be a single .itp file besides forcefield.itp "
        err += f"in {gromacs_param_dir}/toppar/\n"
        raise ValueError(err)
    molec_itp_file = itp_files[0]

    # Look for the forcefield itp file
    ff_itp_file = os.path.join(gromacs_param_dir,"toppar","forcefield.itp")
    if not os.path.exists(ff_itp_file):
        err = "There should be a forcefield.itp file in "
        err += f"in {gromacs_param_dir}/toppar/\n"
        raise ValueError(err)


    # -------------------------------------------------------------------------
    # Create the output directory

    counter = 0
    output_directory = os.path.join(working_dir,
                                    f"{prefix}_custom-param-dir_{counter}")
    while os.path.exists(output_directory):
        counter += 1
        output_directory = os.path.join(working_dir,
                                        f"{prefix}_custom-param-dir_{counter}")

    os.mkdir(output_directory)

    # Copy .gro and forcefield.itp file into the output directory
    shutil.copy(gro_file,os.path.join(output_directory,f"molecule.gro"))
    shutil.copy(ff_itp_file,os.path.join(output_directory,f"forcefield.itp"))

    # -------------------------------------------------------------------------
    # Sync resid names between the .gro and  molecule .itp file. This is required 
    # because the .gro residue column is shorter than the allowed residue name
    # length in the .itp file. Clean up by whacking down the residue names in
    # the .itp file to match what they are called in the .gro file. 

    # Look for unique names in the .gro file
    f = open(gro_file)
    lines = f.readlines()[2:-1]
    f.close()
    unique_resid_in_gro = set([l[5:10] for l in lines])

    # Look for unique residue names in the .itp file
    itp_atoms = ITPManager(molec_itp_file).fields["atoms"]
    unique_resid_in_itp = []
    for atom in itp_atoms:
        unique_resid_in_itp.append(atom[0])
    unique_resid_in_itp = set(unique_resid_in_itp)

    # Find residue names in .gro file that are not in the .itp. 
    gro_mismatches = unique_resid_in_gro.difference(unique_resid_in_itp)
    itp_mismatches = unique_resid_in_itp.difference(unique_resid_in_gro)
    to_gro_name = {}
    for g in gro_mismatches:
        for p in itp_mismatches:

            # Grab the first len(g) letters of the itp name and see if they 
            # match the .gro name. For example "ECLIPA" from the .itp file 
            # would be "ECLIP" in the .gro file. This would mean that 
            # "ECLIPA"[:len("ECLIP")] == "ECLIP". If we had both ECLIPA and 
            # ECLIPB in the .itp file, we would be doomed because both map to 
            # ECLIP...
            if p[:len(g)] == g:
                try:
                    to_gro_name[p]
                    err = f"itp residue {p} ambiguous. maps to {to_gro_name[p]} and {g}\n"
                    raise ValueError(err)
                except KeyError:
                    to_gro_name[p] = g + " "*(len(p) - len(g))

    # Go through the ITP file and rename any atom with mismatched residue names
    # identified above. 
    out_lines = []
    in_atoms = False
    with open(molec_itp_file,'r') as f:

        for line in f:

            if line.startswith("[ atoms ]"):
                in_atoms = True
                out_lines.append(line)
                continue

            if in_atoms:
                if line.startswith("["):
                    in_atoms = False
                else:
                    for p in to_gro_name:
                        line = re.sub(p,to_gro_name[p],line)

            out_lines.append(line)

    # Write molecule .itp file in output directory 
    molec_itp = os.path.join(output_directory,"molecule.itp")
    f = open(molec_itp,'w')
    f.write("".join(out_lines))
    f.close()

    # -------------------------------------------------------------------------
    # extract any positional restraint macros defined in the mdp files. 

    # Get all lines with -DPOSRES from the mdp files
    mdp_posres = _look_for_pattern_in_files("-DPOSRES",
                                             file_search_pattern="*.mdp",
                                             search_dir=gromacs_param_dir)

    all_macro_defs = []
    for hit in mdp_posres:
        line_with_hit = hit[-1].string
        macro_defs = "=".join(line_with_hit.split("=")[1:])
        for m in macro_defs.split():
            all_macro_defs.append(m.strip())

    all_macro_defs = list(set(all_macro_defs))
    all_macro_defs.remove("-DPOSRES")
    
    with open(os.path.join(output_directory,"posres.txt"),"w") as f:
        for m in all_macro_defs:
            f.write(f"{m}\n")

    return output_directory

def merge_ff_itp(imported_param_dir,
                 working_dir):
    """
    Take the forecefield.itp file from a set of custom gromacs parameters and
    merge it into an existing topology file, stripping all duplicate definitions
    from the incoming ff_itp. This is useful because something like charmm-gui
    will generate an itp file for a molecule that includes all necessary 
    atomtypes, etc to describe the molecule.  If we just dropped the .itp file
    into the standard forcefield, we would have duplicate definitions for many 
    parameters. 

    Parameters
    ---------- 
    ff_itp : str
        forcefield itp file with custom parameters
    topo_file : str
        current system topology file (which likely brings in forcefield itp
        files via #include directories)

    Returns
    -------
    itp_file : str
        newly written out itp file
    """
    
    ff_itp = os.path.join(imported_param_dir,"forcefield.itp")
    topo_file = os.path.join(working_dir,"topol.top")

    # Load custom itp file into an ITPManager
    custom_param_ff = ITPManager(ff_itp)

    # Get all itp files referenced by the topo file
    itp_files = _get_included_files(some_itp=topo_file,
                                    base_path=working_dir)

    # Mask any entries from the new itp file that duplicate entries already 
    # present in the topology
    for f in itp_files:
        custom_param_ff.mask_existing_entries(other_itp_file=f)
    
    # Mask defaults -- this will always conflict and can't be imported after 
    # forcefield #include directive without error
    custom_param_ff.mask_field(field="defaults")

    # Find a new itp file name
    itp_counter = 0
    custom_itp_file = f"custom-param-ff_{itp_counter}.itp"
    while os.path.exists(os.path.join(working_dir,custom_itp_file)):
        itp_counter += 1
        custom_itp_file = f"custom-param-ff_{itp_counter}.itp"

    # Write to the new itp file name
    itp_file_to_write = os.path.join(working_dir,custom_itp_file)
    custom_param_ff.write_itp(itp_file_to_write)
    
    # Add new #include statement pointing to new itp file
    _insert_new_include(what_to_include=custom_itp_file,
                        itp_file=topo_file,
                        insert_index=0,
                        insert_before=False)

    # Return the name of the file we wrote out
    return itp_file_to_write

def merge_molec_itp(imported_param_dir,
                    working_dir,
                    ignore_posres=False):
    """
    Take the .itp file describing the molecule atoms etc. and add it into the
    working directory. Check to see if the molecule.itp is already included in 
    the topology. If it is not, import the molecule and add an #include 
    statement to the topology. Then add another instance of the molecule to 
    the topol.top [ molecules ] list. 

    Parameters
    ----------
    imported_param_dir : str
        directory output from import_custom_params
    working_dir : str
        current working directory
    ignore_posres : bool, default=False
        ignore any positional restraint information
    """

    molec_itp = os.path.join(imported_param_dir,"molecule.itp")
    topol_top = os.path.join(working_dir,"topol.top")

    im = ITPManager(molec_itp)
    moleculetype = list(im.fields["moleculetype"].keys())[0]

    # See if this molecule type has already been defined somewhere
    need_to_include_itp = True
    included_itp_files = _get_included_files(topol_top)
    for other_itp in included_itp_files:
        other_im = ITPManager(other_itp)
        if "moleculetype" in other_im.fields:
            other_moleculetype = list(other_im.fields["moleculetype"].keys())[0]
            if other_moleculetype == moleculetype:
                need_to_include_itp = False
                break

    if need_to_include_itp:
        
        # Find a new itp file name
        itp_counter = 0
        custom_itp_file = f"custom-param-molec_{itp_counter}.itp"
        while os.path.exists(os.path.join(working_dir,custom_itp_file)):
            itp_counter += 1
            custom_itp_file = f"custom-param-molec_{itp_counter}.itp"

        shutil.copy(molec_itp,os.path.join(working_dir,custom_itp_file))

        # If we're ignoring position restraints, remove them from the itp file
        if ignore_posres:

            # Load itp to modify into an ITPManager
            itp_to_mod = os.path.join(working_dir,custom_itp_file)
            im = ITPManager(itp_to_mod)
            
            # Mask any field ending with _restraints
            for field in im.fields:
                cols = field.split("_")
                if cols[-1] == "restraints":
                    im.mask_field(field)
            
            # Write masked itp
            im.write_itp(itp_to_mod)

        # Add new #include statement pointing to new itp file
        _insert_new_include(what_to_include=custom_itp_file,
                            itp_file=topol_top,
                            insert_index=-1,
                            insert_before=False)

    # Append the molecule to the end of the [ molecules ] field in the 
    # topol.top file
    with open(topol_top) as f:
        topol_lines = f.readlines()

    in_molecules = False
    for i, line in enumerate(topol_lines):
        if line.startswith("[ molecules ]"):
            in_molecules = True
            continue

        if in_molecules and line.startswith("["):
            break

    topol_lines.insert(i+1,f"{moleculetype}    1\n")
    
    with open(topol_top,"w") as f:
        f.writelines(topol_lines)    
       
def _pdb_into_gro(pdb_file,
                  gro_file):
    """
    Load the coordinates from a pose pdb file into a gro file.  
   
    Parameters
    ---------- 
    pdb_file : str
        pdb file with new atom coordinates.  these atoms must be in the same
        order and have the same first four characters in the names as the 
        custom parameter gro file.  This is what is written out by pymol if you
        load parameter .gro file and then export with the 'Original atom order'
        and 'Retain atom ids' options checked. 
    gro_file : str
        gro file that will have data loaded into it.

    Returns
    -------
    out_lines : list
        list of lines of the new gro file
    """
    
    # Read pdb file, grabbing atom names and coordinates in nm (note conversion
    # from angstroms!)
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line[:6] in ["ATOM  ","HETATM"]:
                atom_name = line[12:16].strip()
                coord = [float(line[30+8*i:(38+8*i)])/10.0 for i in range(3)]
                atoms.append((atom_name,coord))
        
    # Go through gro file and replace atom coordinates with those from the pdb
    # file
    out_lines = []
    with open(gro_file) as f:
        for i, line in enumerate(f):
            
            # Skip two header lines
            if i < 2:
                out_lines.append(line)
                continue
            
            # Last line is not atom in gro file
            if i - 2 == len(atoms):
                out_lines.append(line)
                continue
            
            # Make sure atoms match
            atom = line[10:15].strip()
            if atom[:4] != atoms[i-2][0]:
                err = "apparent mismatch in atoms in gro and pdb file\n"
                err += f"   {atom[:4]} vs {atoms[i-2][0]}\n"
                raise ValueError(err)
            
            # Load in new coordinates
            front = line[:20]
            coord = "".join(["{:>8.3f}".format(c) for c in atoms[i-2][1]])
            back = line[44:]
            
            out_lines.append(f"{front}{coord}{back}")
            
    return out_lines

def merge_molec_coord(coord_file,
                      imported_param_dir,
                      working_dir):
    """
    Append a set of molecular coordinates to system.gro. 

    Parameters
    ----------
    coord_file : str
        file with the coordinates of the molecule
    imported_param_dir : str
        directory output from import_custom_params that corresponds to the 
        molecule in coord_file
    working_dir : str
        current working directory
    """

    if coord_file[-4:].lower() == ".pdb":
        new_gro_lines = _pdb_into_gro(pdb_file=coord_file,
                                      gro_file=os.path.join(imported_param_dir,"molecule.gro"))
    elif coord_file[-4:].lower() == ".gro":
        with open(coord_file) as f:
            new_gro_lines = f.readlines()
    else:
        err = "coord_file should have .pdb or .gro extension\n"
        raise ValueError(err)

    new_gro_lines = new_gro_lines[2:-1]

    system_gro = os.path.join(working_dir,"system.gro")
    with open(system_gro) as f:
        gro_lines = f.readlines()

    out_lines = gro_lines[2:-1]

    # Figure out residue offset we need to add to new gro file
    this_start = int(new_gro_lines[0][:5])
    if len(out_lines) == 0:
        last_resid = 0
    else:
        last_resid = int(out_lines[-1][:5])
    offset = last_resid - this_start + 1

    # Go through lines from new gro file and stick into system gro file
    for g in new_gro_lines:
        v = int(g[:5]) + offset
        out_lines.append("{:5d}{}".format(v,g[5:]))

    # Construct finalized gro file
    final_gro_contents = [gro_lines[0]]
    final_gro_contents.append("{:d}\n".format(len(out_lines)))
    final_gro_contents.extend(out_lines)
    final_gro_contents.append(gro_lines[-1])

    f = open(system_gro,'w')
    f.write("".join(final_gro_contents))
    f.close()

def merge_posres(imported_param_dir,
                 working_dir):
    """
    Take any -DPOSRES from the imported parameters and append them to the 
    positional restraint calls in any mdp files. 

    Parameters
    ----------
    imported_param_dir : str
        directory output from import_custom_params
    working_dir : str
        current working directory
    """
    
    with open(os.path.join(imported_param_dir,"posres.txt")) as f:
        posres = [p.strip() for p in f.readlines()]

    posres_lines = _look_for_pattern_in_files(text_search_pattern="-DPOSRES",
                                              file_search_pattern="*.mdp",
                                              search_dir=working_dir)
    
    for hit in posres_lines:
        
        filename = hit[0]
        line_number = hit[1]

        with open(filename) as f:
            mdp_lines = f.readlines()

        line = mdp_lines[line_number]

        # Deal with possible comments on the line
        split_on_comment = line.split(";")
        pre_comment = split_on_comment[0]
        if len(split_on_comment) > 1:
            post_comment = " ;" + ";".join(split_on_comment[1:])
        else:
            post_comment = "\n"

        existing_defs = "=".join(pre_comment.split("=")[1:])
        
        new_line = ["define","=",existing_defs]
        new_line.extend(posres)
        new_line.append(post_comment)

        mdp_lines[line_number] = " ".join(new_line)

        with open(filename,"w") as f:
            f.writelines(mdp_lines)

def setup_md(output_dir,
             standard_pdb_file=None,
             custom_poses=None,
             custom_params=None,
             ff_source_dir=None,
             run_files_dir=None,
             overwrite=False,
             box_size=1.0,
             ion_conc=0.1,
             pos_charge=1,
             neg_charge=-1,
             pos_ion="NA",
             neg_ion="CL",
             water_model="tip3p",
             use_resid_offset=True,
             ignore_posres=False):
    """
    Set up a GROMACS MD run, possibly with custom forcefield parameters
    generated by charmm-gui. This will add a solvent box and ions to neutralize
    the system.

    Parameters
    ----------
    output_dir : str
        directory to write files to
    ff_source_dir : str
        directory of forcefield in gromacs format (e.g. CHARMM36m.ff)
    run_files_dir : str
        directory with mdp and other input files used to do a run

    standard_pdb_file : str, optional
        pdb file containing all of the standard molecules you wish to include
        in the simulation. this corresponds to protein, nucleic acids, many ions,
        and some lipids & carbohydrates. Basically, anything that can be 
        handled automatically with your ff_source_dir and pdb2gmx. 
    custom_poses : str, optional
        pdb or gro file (or list of such files) with poses for molecules that 
        require custom parameters. 
    custom_params : str
        gromacs dir (or list of such directories) spit out by charmm-gui
        corresponding to the coordinates in custom_poses

    overwrite : bool, default=False
        overwrite output_dir

    box_size : float, default=1.0
        length of edge of box away from molecule in all directions (nm)
    ion_conc : float, default=0.1
        ion concentration in M.  will ions to this conc in numbers that will
        neutralize system. 
    pos_charge : int, default=1
        charge of the positive ions to add
    neg_charge : int, default=-1
        charge of the negative ions to add
    pos_ion : str, default="NA"
        name of the positive ions to add
    neg_ion : str, default="CL"
        name of the negative ions to add
    water_model : str, default="tip3p"
        water model to use. should be one of: spc, spce, tip3p, tip4p, tip5p

    use_resid_offset : bool, default=True
        offset chains by 1000 residues for any chains in standard_pdb_file. 
        chain A would be 1-100, chain B 1001-1100, etc. This makes it easier to
        keep track of chains after the simulation.
    ignore_posres : bool, default=False
        ignore positional restraints coming in via the custom parameters
              
    WARNING: this does *not* check for forcefield consistency.  You could 
    generate the custom ff on charmm-gui using CHARMM36m and then stick it into
    the AMBER03 forcefield.  This script may work, without error, but could 
    lead to wacky outcomes in the simulations. 
    """

    # -------------------------------------------------------------------------
    # Check arguments 

    # Create a list of custom pdb files
    if custom_poses is None:
        custom_poses = []
    if issubclass(type(custom_poses),str):
        custom_poses = [custom_poses]

    # Create a list of custom parameters
    if custom_params is None:
        custom_params = []
    if issubclass(type(custom_params),str):
        custom_params = [custom_params]

    # Make sure user specified at least one standard or custom pdb file
    if len(custom_poses) == 0 and standard_pdb_file is None:
        err = "You must specify at least one standard or custom pdb file.\n"
        raise ValueError(err)
    
    # Make sure there are the same number of custom pdb files and parameters
    if len(custom_poses) != len(custom_params):
        err = "Each custom pdb file must have an associated custom param directory\n"
        raise ValueError(err)

    # Make sure all standard pdb files are files
    if standard_pdb_file is not None:
        if not os.path.isfile(standard_pdb_file):
            err = f"standard_pdb_file {standard_pdb_file} is not a file\n"
            raise FileNotFoundError(err)

    # Make sure all custom poses are files
    if len(custom_poses) > 0:
        for f in custom_poses:
            if not os.path.isfile(f):
                err = f"custom_pdb_file {f} is not a file\n"
                raise FileNotFoundError(err)     

    # Make sure all custom parameters are directories
    if len(custom_params) > 0:
        for p in custom_params:
            if not os.path.isdir(p):
                err = f"custom_param {p} is not a directory\n"
                raise FileNotFoundError(err)

    # Check to see if output directory exists and whether to overwrite
    if os.path.isdir(output_dir):
        if overwrite:
            shutil.rmtree(output_dir)
        else:        
            err = f"'{output_dir}' already exists!\n"
            raise FileExistsError(err)


    if ff_source_dir is None:

        script_dir = os.path.dirname(os.path.realpath(__file__))
        ff_source_dir = os.path.join(script_dir,DEFAULT_FF_DIR)
        if not os.path.isdir(ff_source_dir):
            err = "\n\nIf no ff_source_dir is specified, setup_md.py looks for the\n"
            err += f"directory '{DEFAULT_FF_DIR}' inside the directory that \n"
            err += "contains setup_md.py. This directory was not found in\n"
            err += f"'{script_dir}'.\n" 
            err += "Either specify the ff_source_dir manually or make sure that\n"
            err += f"'{DEFAULT_FF_DIR}' is in the same directory\n"
            err += "as the setup_md.py script.\n\n"
            raise FileNotFoundError(err)

    if not os.path.isdir(ff_source_dir):
        err = f"ff_source_dir {ff_source_dir} is not a directory\n"
        raise ValueError(err)

    if run_files_dir is None:

        script_dir = os.path.dirname(os.path.realpath(__file__))
        run_files_dir = os.path.join(script_dir,DEFAULT_RUN_FILES_DIR)
        if not os.path.isdir(run_files_dir):
            err = "\n\nIf no run_files_dir is specified, setup_md.py looks for the\n"
            err += f"directory '{DEFAULT_RUN_FILES_DIR}' inside the directory that \n"
            err += "contains setup_md.py. This directory was not found in\n"
            err += f"'{script_dir}'.\n"
            err += "Either specify the run_files_dir manually or make sure that\n"
            err += f"'{DEFAULT_RUN_FILES_DIR}' is in the same directory\n"
            err += "as the setup_md.py script.\n\n"
            raise FileNotFoundError(err)

    if not os.path.isdir(run_files_dir):
        err = f"run_files_dir {run_files_dir} is not a directory\n"
        raise ValueError(err)


    # -------------------------------------------------------------------------
    # Populate the source and working directories with with user inputs

    # Create output directory 
    os.mkdir(output_dir)
    
    # Create a source directory to store the raw ff and charmm-gui output
    src_dir = os.path.join(output_dir,"src")
    os.mkdir(src_dir)
    
    # Create working directory
    working_dir = os.path.join(output_dir,"working")
    os.mkdir(working_dir)
    
    # Create final direcotry
    final_dir = os.path.join(output_dir,"final")
    os.mkdir(final_dir)
    
    # Copy in the forcefield directory
    ff_source_dir_base = os.path.basename(ff_source_dir)
    shutil.copytree(ff_source_dir,os.path.join(src_dir,ff_source_dir_base))
    shutil.copytree(ff_source_dir,os.path.join(working_dir,ff_source_dir_base))
    
    # Copy in the mdp files
    run_files_dir_base = os.path.basename(run_files_dir)
    shutil.copytree(run_files_dir,os.path.join(src_dir,run_files_dir_base))
    run_files = glob.glob(os.path.join(run_files_dir,"*.*"))
    for f in run_files:
        shutil.copy(f,working_dir)
    
    tmp_base = "temporary-setup-md-file"
    
    # -------------------------------------------------------------------------
    # Create system topology and gro files

    print("Creating initial topology",flush=True)    
    initialize_directory(standard_pdb_file=standard_pdb_file,
                         working_dir=working_dir,
                         water_model=water_model,
                         ff_source_dir=ff_source_dir_base,
                         use_resid_offset=use_resid_offset,
                         prefix=f"{tmp_base}_initialize")
    
    # -------------------------------------------------------------------------
    # Bring in any custom molecules

    for i in range(len(custom_params)):

        prefix = f"{tmp_base}_custom-{i}"

        # Get the parameter dir
        param_dir = custom_params[i]

        # Copy the custom parameters directory into the source
        gromacs_param_dir_base = os.path.basename(os.path.dirname(param_dir))
        gromacs_param_dir = os.path.join(src_dir,
                                         f"custom_params_{i}",
                                         gromacs_param_dir_base)
        shutil.copytree(param_dir,gromacs_param_dir)
        
        # Import custom parameter file
        imported_param_dir = import_custom_params(gromacs_param_dir=gromacs_param_dir,
                                                  working_dir=working_dir,
                                                  prefix=prefix)
        
        # Merge the forcefield.itp file with the existing .itp files in the
        # directory
        merge_ff_itp(imported_param_dir=imported_param_dir,
                     working_dir=working_dir)

        # Add the molecule to the topol.top file
        merge_molec_itp(imported_param_dir=imported_param_dir,
                        working_dir=working_dir,
                        ignore_posres=ignore_posres)
        
        # Add the molecule coordinates to the system.gro file
        pose_name = f"custom_pose_{i}.{custom_poses[i][-3:]}"
        pose_file = os.path.join(src_dir,pose_name)
        shutil.copy(os.path.abspath(custom_poses[i]),pose_file)
        merge_molec_coord(coord_file=pose_file,
                          imported_param_dir=imported_param_dir,
                          working_dir=working_dir)
        
        # Merge the position restraints detected into any .mdp files unless we
        # are ignoring position restraints
        if not ignore_posres:
            merge_posres(imported_param_dir=imported_param_dir,
                        working_dir=working_dir)

    # Create an .ndx file that has the user inputs, not solvent ions or waters
    # added later
    create_pre_solvate_ndx(working_dir=working_dir,
                           prefix=f"{tmp_base}_pre-solvate-idx")

    # -------------------------------------------------------------------------
    # Finalize the system by adding a solvent box, waters, and counter ions to 
    # neutralize

    solvate_and_ionize(working_dir=working_dir,
                       tmp_base=tmp_base,
                       box_size=box_size,
                       ion_conc=ion_conc,
                       pcharge=pos_charge,
                       ncharge=neg_charge,
                       pname=pos_ion,
                       nname=neg_ion)

    # Get list of non-temporary files
    non_tmp = [c for c in os.listdir(working_dir) if not c.startswith(tmp_base)]
    non_tmp = [c for c in non_tmp if not c.startswith("#")]
    
    # Remove a couple of temporary files that were not caught by above
    to_remove = ["mdout.mdp","ions.tpr","ions.mdp"]
    non_tmp = [c for c in non_tmp if c not in to_remove]

    # Copy files into the output directory
    for c in non_tmp:        
        source = os.path.join(working_dir,c)
        target = os.path.join(final_dir,c)
        if os.path.isdir(source):
            shutil.copytree(source,target)
        else:
            shutil.copy(source,target)

def main(argv=None):
    """
    Run the setup_md function from the command line. 

    Parameters
    ----------
    argv : list, optional
        list of arguments. if not specified, pull sys.argv[1:]
    """

    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(prog='setup_md.py',
                                     description='prepare a GROMACS simulation, possibly with custom parameters generated by charmm-gui')
    
    parser.add_argument("output_dir",
                        help="Directory to write files to. The directory in which you should run the simulation will be output_dir/final")
    
    parser.add_argument("--overwrite",
                        action="store_true",
                        help="overwrite existing output directory if it exists")
    parser.add_argument("--version",
                        action="version",
                        version=f"%(prog)s {__version__}")
    
    input_group = parser.add_argument_group("Molecular Inputs",description="Specify what molecules are going to be used in the simulation.")

    input_group.add_argument("--standard",
                             nargs=1,
                             type=str,
                             help="Single pdb file containing all standard molecules in the simulation (proteins, nucleic acids, standard ions, etc.)")
    input_group.add_argument("--custom",
                             nargs="+",
                             help="pdb or gro file(s) with custom molecules. If these are specified, a --params argument must be specified for each molecule.")
    input_group.add_argument("--params",
                             nargs="+",
                             type=str,
                             help="Location of 'gromacs' directory(s) generated by charmm-gui. One must be specified for each --custom pdb")

    sim_group = parser.add_argument_group("Forcefield & Run Inputs",description="Specify forcefield and simulation input files for the simulation.")
    sim_group.add_argument("--ff-source-dir",
                           nargs=1,
                           type=str,
                           help=f"Directory of forcefield in gromacs format. Default is to look in the directory containing setup_md.py for '{DEFAULT_FF_DIR}'")
    sim_group.add_argument("--run-files-dir",
                           nargs=1,
                           type=str,
                           help=f"Directory with mdp and other files needed files to set up run. Default is to look in the directory containing setup_md.py for '{DEFAULT_RUN_FILES_DIR}'")
    
    solv_group = parser.add_argument_group("Solvation Inputs",description="Options controlling how the molecules will be solvated.")
    solv_group.add_argument("--box-size",
                            default=1.0,
                            type=float,
                            help="Length of edge of box away from molecule in all directions (nm). Default is 1.0")
    solv_group.add_argument("--ion-conc",
                            default=0.1,
                            type=float,
                            help="Ion concentration in M. Script will add ions to this conc in numbers to neutralize system. Default is 0.1.")
    solv_group.add_argument("--pos-charge",
                            default=1,
                            type=int,
                            help="Charge of the positive ions to add. Default is 1")
    solv_group.add_argument("--neg-charge",
                            default=-1,
                            type=int,
                            help="Charge of the negative ions to add. Default is -1")
    solv_group.add_argument("--pos-ion",
                            default="NA",
                            type=str,
                            help="Name of the positive ions to add. Default is 'NA'")
    solv_group.add_argument("--neg-ion",
                            default="CL",
                            type=str,
                            help="Name of the negative ions to add. Default is 'CL'")
    solv_group.add_argument("--water-model",
                            default="tip3p",
                            help="Water model to use. Should be one of: spc, spce, tip3p, tip4p, tip5p. Default is 'tip3p' ")

    
    adv_control = parser.add_argument_group("Advanced Control",description="Options for controlling how script treats inputs.")
    adv_control.add_argument("--no-resid-offset",
                             action='store_true',
                             help="By default, this script will offset each chain in the standard pdb file by 1000 residues (e.g. A will be 1-100, B 1001-1100). This flag turns this behavior off.")
    adv_control.add_argument("--ignore-posres",
                             action="store_true",
                             help="Ignore any positional restraint data in the custom parameter files. The GROMACS positional restraints coming out of charmm-gui are sometimes formatted incorrectly. This causes GROMACS to crash during a position-restrained run, complaining about incorrectly defined constraints. This fixes the problem by turning off the constraints for custom molecules. (This should be fine for small molecule ligands, particularly bound in a protein pocket that will have constraints on during equilibration).")

    parsed = parser.parse_args()

    # Build keyword arguments
    kwargs = {}
    kwargs["output_dir"] = parsed.output_dir
    kwargs["ff_source_dir"] = parsed.ff_source_dir
    kwargs["run_files_dir"] = parsed.run_files_dir
    
    if parsed.standard is None:
        kwargs["standard_pdb_file"] = None
    else:
        kwargs["standard_pdb_file"] = parsed.standard[0]
    kwargs["custom_poses"] = parsed.custom
    kwargs["custom_params"] = parsed.params

    kwargs["overwrite"] = parsed.overwrite
    kwargs["box_size"] = parsed.box_size
    kwargs["ion_conc"] = parsed.ion_conc
    kwargs["pos_charge"] = parsed.pos_charge
    kwargs["neg_charge"] = parsed.neg_charge
    kwargs["pos_ion"] = parsed.pos_ion
    kwargs["neg_ion"] = parsed.neg_ion
    kwargs["water_model"] = parsed.water_model

    if parsed.no_resid_offset:
        kwargs["use_resid_offset"] = False
    else:
        kwargs["use_resid_offset"] = True

    kwargs["ignore_posres"] = parsed.ignore_posres
    

    setup_md(**kwargs)


if __name__ == "__main__":
    main()  
