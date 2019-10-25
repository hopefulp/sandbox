#!/opt/applic/epd/bin/python

# python modules
import sys
import os
import string
import copy
import math
import re
import time

version = '20160523'

def getLogData(log_file, out_file, requested_key, silent=True):
    """
getLogData: GET thermodynamic data from a LAMMPS log file
    if out_file is specified, the output will be written on the file.
    or, the list of [step value1 value2 ...] will be returned.
    """


    # initialization
    n_count = 0; TICK = 10000;


    # parse section and extract all thermodynamic data
    pat1 = re.compile(r"Step\s*(\S*)\s")
    pat2 = re.compile(r"([A-z]+)\s+=\s*-*\d+\.\d+")    # keywords
    pat3 = re.compile(r"[A-z]+\s+=\s*(-*\d+\.\d+)")    # numbers


    # log file loading
    myLOG = open(log_file)
    myLOG_output = [];


    # preprocessing
    if not silent: print("\n== Step 1. Preprocessing the LAMMPS log file")
    wc_log_file = os.popen("wc -l " + log_file).read()    # file newline counts
    log_file_length = int(wc_log_file.split()[0]); str_log_file_length = str(log_file_length);
    t1 = time.time(); t2 = 0;
    while 1:
        try:
            line = myLOG.readline()

            if not line:
                break;

            line = line.replace("\n", "")

            n_count += 1;
            if not silent:
                if n_count % TICK == 0 or n_count == log_file_length:
                    # estimated time calculation
                    t2 = time.time()
                    elapsed = t2 - t1
                    estimated = elapsed * ((log_file_length - n_count) // TICK)    

                    # string conversion
                    str_n_count = str(n_count)
                    str_estimated = "{0:4.1f}".format(estimated)
                    sys.stdout.write("\rPreprocessing.. "+ str_n_count + " / " + str_log_file_length + ".. " + "(" + str_estimated + " sec remain)")
                    sys.stdout.flush()

                    # reassigning time
                    t1 = t2;

            # keyword checking
            """    ### ORIGINAL version
            if "angle_coeff" in line or "angle_style" in line or "atom_modify" in line or "atom_style" in line or \
                "bond_coeff" in line or "bond_style" in line or "boundary" in line or "change_box" in line or \
                "variable" in line or "data" in line or "group" in line or "dump" in line or \
                "clear" in line or "communicate" in line or "compute" in line or "compute_modify" in line or \
                "create_atoms" in line or "create_box" in line or "delete_atoms" in line or "delete_bonds" in line or \
                "dielectric" in line or "dihedral_coeff" in line or "dihedral_style" in line or "dimension" in line or \
                "displace_atoms" in line or "displace_box" in line or "dump" in line or "dump" in line or "image" in line or \
                "dump_modify" in line or "echo" in line or "fix" in line or "fix_modify" in line or \
                "group" in line or "if" in line or "improper_coeff" in line or "improper_style" in line or \
                "include" in line or "jump" in line or "kspace_modify" in line or "kspace_style" in line or \
                "label" in line or "lattice" in line or "log" in line or "mass" in line or \
                "minimize" in line or "min_modify" in line or "min_style" in line or "neb" in line or \
                "neigh_modify" in line or "neighbor" in line or "newton" in line or "next" in line or \
                "package" in line or "pair_coeff" in line or "pair_modify" in line or "pair_style" in line or \
                "pair_write" in line or "partition" in line or "prd" in line or "print" in line or \
                "processors" in line or "read_data" in line or "read_restart" in line or "region" in line or \
                "replicate" in line or "reset_timestep" in line or "restart" in line or "run" in line or \
                "run_style" in line or "set" in line or "shell" in line or "special_bonds" in line or \
                "suffix" in line or "tad" in line or "temper" in line or "thermo" in line or \
                "thermo_modify" in line or "thermo_style" in line or "timestep" in line or "uncompute" in line or \
                "undump" in line or "unfix" in line or "units" in line or "variable" in line or \
                "velocity" in line or "write_restart" in line:
                pass;
            """

            ### 2nd version
            if "coeff" in line or "modify" in line or "style" in line or \
                "boundary" in line or "change_box" in line or \
                "variable" in line or "data" in line or "group" in line or "dump" in line or \
                "clear" in line or "communicate" in line or "compute" in line or \
                "create" in line or "delete" in line or \
                "dielectric" in line or "dimension" in line or \
                "displace_atoms" in line or "displace_box" in line or "image" in line or \
                "echo" in line or "fix" in line or \
                "group" in line or "if" in line or \
                "include" in line or "jump" in line or \
                "label" in line or "lattice" in line or "log" in line or "mass" in line or \
                "minimize" in line or "neb" in line or \
                "neighbor" in line or "newton" in line or "next" in line or \
                "package" in line or \
                "write" in line or "partition" in line or "prd" in line or "print" in line or \
                "processors" in line or "read" in line or "region" in line or \
                "replicate" in line or "time" in line or "restart" in line or "run" in line or \
                "set" in line or "shell" in line or "special_bonds" in line or \
                "tad" in line or "temper" in line or "thermo" in line or \
                "uncompute" in line or \
                "units" in line or "variable" in line or \
                "velocity" in line:
                pass;

            elif "Step" in line or "TotEng" in line or "PotEng" in line or "E_dihed" in line or "E_coul" in line or \
                "Temp" in line or "Volume" in line:
                myLOG_output.append(line)

            elif "KinEng" in line or "E_bond" in line or "E_impro" in line or "E_long" in line or \
                "E_angle" in line or "E_vdwl" in line or "Press" in line or  \
                "v_" in line or "c_" in line:
                myLOG_output.append(line)

        except KeyboardInterrupt:
            nu.die("Keyboard Break.. Exiting.")


    # find the number of last step
    templist = copy.deepcopy(myLOG_output)
    laststep = "";
    while 1:
        l = templist.pop()
        if "Step" in l:
            laststep = l
            break;
    laststep = int(re.search(pat1, laststep).group(1))    # last step
    str_laststep = str(laststep)

    if not silent: sys.stdout.write("\nDone.\n")


    if not silent: print("\n== Step 2. Getting Thermodynamic Properties")

    section = ""; thermo_data = []; line = ""; lines = [];

    # first line: pass!!
    section = myLOG_output[0]
    myLOG_output = myLOG_output[1:]

    t1 = time.time(); t2 = 0;


    # parse and check requested keywords
    if not requested_key == "":
        requested_key = requested_key.split();
    

    # treat preprocessed log data in myLOG_output
    try:
        for line in myLOG_output:
    
            section_full = False;
    
            if "Step" in line:
                lines.append(section)
                section = line;
                section_full = False;
            else:
                section += line;

        lines.append(section)

        # actual parsing
        for section in lines:
            step = int(re.search(pat1, section).group(1)) # steps
            keyword = re.findall(pat2, section) # keywords
            value = re.findall(pat3, section) # numbers

            # display step progress
            if not silent:
                if step % TICK == 0 or step == laststep:
                    # estimated time calculation
                    t2 = time.time()
                    elapsed = t2 - t1
                    estimated = elapsed * ((laststep - step) // TICK)    

                    str_step = str(step)
                    str_estimated = "{0:4.1f}".format(estimated)
                    sys.stdout.write("\rParsing the Step " + str_step + " / " + str_laststep + ".. " + "(" + str_estimated + " sec remain)")
                    sys.stdout.flush()

                    # reassign time
                    t1 = t2;

            if len(keyword) != len(value):
                nu.die("Suspicious parsing at " + str(step) + ".. Exiting.")

            # get values
            temp = [];
            if requested_key == "":
                # no specified keywords: return the whole data [key1 value1 key2 value2 ...]
                temp.append("Step"); 
                temp.append(str(step))

                for i in range(0, len(keyword)):
                    temp.append(str(keyword[i]));
                    temp.append(str(value[i]))
            else:
                temp.append(str(step))
                # specified keywords
                for rk in requested_key:
                    try:
                        temp.append(value[keyword.index(rk)])
                    except:
                        temp.append("NA")

            thermo_data.append(temp)

    except KeyboardInterrupt:
        nu.die("Keyboard Break.. Exiting.")

    if not requested_key == "":
        thermo_data.insert(0, ["Step"] + requested_key)


    # if out_file is specified, then the result will be written in the file.
    if out_file != "":
        f_out_file = open(out_file, "w")        # output file
        for i in thermo_data:
            output = "\t".join(i);
            output += "\n"
            f_out_file.write(output)
        if not silent: print("\nParsed thermodynamic data is saved on " + out_file + " ..")
        f_out_file.close()

    if not silent: sys.stdout.write("Done.\n")

    return thermo_data;    # success

    ### end of getLogData


def getTrjInfo(trj_file, silent=True):

    if not silent: print("Looking up trajectory file..")
    
    def get_line(trj_file):
        """
        Open a lazy file stream. This is considered to increase a file reading speed.
        """
        with open(trj_file) as file:
            for i in file:
                yield i

    # Initialize
    timesteps = []; 

    # ToDo:
    # Find a hidden info file. filename: .$trjfile.info
    # if modtime of trj and modtime in infofile is same:
    #       return the timestep in infofile
    # or proceed

    def scan():
        header = [];
        for line in get_line(trj_file):
            if "ITEM: ATOMS" in line:
                header.append(line)
                return header
            else:
                header.append(line)

    n_header = len(scan())
    n_atoms = int(scan()[3])

    # get timesteps
    if not silent: print("Reading LAMMPS trajectory " + str(trj_file) + " and scanning timesteps.. this may take some time.")

    dumpatom = get_line(trj_file)

    while 1:
        try:
            chunk = [next(dumpatom) for i in range(n_atoms+n_header)]
        except StopIteration:
            break;

        timestep = int(chunk[1])
        timesteps.append(timestep)

        sys.stdout.write('\r' + "Fetched timestep: " + str(timestep))
        sys.stdout.flush()

    # write information to infofile
    # in the file: 1) modification time of the trajectory file, 2) n_timestep

    sys.stdout.write('\n')

    return timesteps, n_header, n_atoms
