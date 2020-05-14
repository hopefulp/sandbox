import re
import os
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

version = '20180327'
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s: %(message)s", datefmt="%y-%m-%_d %H:%M")
logger = logging.getLogger("Log_parser")

pat1 = re.compile(r"Step[-= ]*(\d*)")  # steps
pat2 = re.compile(r"(\S+)\s+=\s*-*\d+\.\d+")  # keywords
pat3 = re.compile(r"\S+\s+=\s*(-*\d+\.\d+)")


def _refine(orig_line):
    orig_line = orig_line.replace("\n", "")
    line = orig_line.lower()

    if "coeff" in line or "modify" in line or "style" in line or \
            "boundary" in line or "change_box" in line or "box" in line or \
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
            "processor" in line or "read" in line or "region" in line or \
            "replicate" in line or "time" in line or "restart" in line or "run" in line or \
            "set" in line or "shell" in line or "special_bonds" in line or \
            "tad" in line or "temper" in line or "thermo" in line or \
            "uncompute" in line or \
            "units" in line or "variable" in line or \
            "velocity" in line or "WARNING" in line or "#" in line:
        return ""

    elif "step " in line or "toteng" in line or "poteng" in line or "e_dihed" in line or "e_coul" in line or \
            "temp" in line or "volume" in line:
        return orig_line

    elif "kineng" in line or "e_bond" in line or "e_impro" in line or "e_long" in line or \
            "e_angle" in line or "e_vdwl" in line or "press" in line or \
            "v_" in line or "c_" in line:
        return orig_line

    elif "density" in line or "lx" in line or "ly" in line or "lz" in line or \
            "cella" in line or "cellb" in line or "cellc" in line or \
            "cellalpha" in line or "cellbeta" in line or "cellgamma" in line or \
            "pxx" in line or "pyy" in line or "pzz" in line or "pxy" in line or "pxz" in line or "pyz" in line or \
            "e_tail" in line or "enthalpy" in line or "e_pair" in line or "e_dihed" in line:
        return orig_line

    else:
        return ""


def lindexsplit(some_list, args):
    args = [0] + args + [len(some_list) + 1]

    # For a little more brevity, here is the list comprehension of the following
    # statements:
    #    return [some_list[start:end] for start, end in zip(args, args[1:])]

    my_list = []
    for start, end in zip(args, args[1:]):
        my_list.append(some_list[start:end])
    return my_list


class Log:
    def __init__(self, log_file, verbose=False):
        assert isinstance(log_file, str)
        if not os.path.exists(log_file):
            raise IOError("LAMMPS log file %s does not exist!" % log_file)

        if verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)

        # preprocessing
        logger.debug("Preprocessing the LAMMPS log file")
        with open(log_file, 'rt') as f:
            reset_timestep = [line_no for line_no, line in enumerate(f.readlines()) if "reset_timestep" in line]

        ans = 0
        if reset_timestep:
            # ans = int(input("'reset_timestep' keyworkd found in the log file. (line %s) \nSelect the part you want to "
            #                 "proceed: (1) before (2) after (3) whole [Default:2] " % reset_timestep) or "2")
            ans = 2

        with open(log_file, 'rt') as f:
            lines = [line for line in f.readlines()]

        if len(reset_timestep) <= 1:
            if ans == 1:
                lines = lines[:reset_timestep[0]]
            elif ans == 2:
                lines = lines[reset_timestep[0]:]
            else:
                pass
        else:
            raise RuntimeError("Too many reset_timestep keyword found!")

        lines = [_refine(line) for line in lines]
        lines = [line for line in lines if line]

        # slice
        step_indexes = []
        for index, line in enumerate(lines):
            if "Step" in line:
                step_indexes.append(index)
        #step_indexes.remove(0)
        print(step_indexes)

        sections = lindexsplit(lines, step_indexes)
        sections = [" ".join(section) for section in sections]

        thermo = {}
        for section in sections:
            try:
                step = int(re.search(pat1, section).group(1))
                keyword = re.findall(pat2, section)
                value = re.findall(pat3, section)
            except ValueError:
                raise ValueError("An error occurred when parsing this line: %s" % section)
            else:
                thermo[step] = dict(zip(keyword, value))

        self.df = pd.DataFrame(thermo, dtype=float).T

    def plot(self, keyword):
        """
        Draw a plot with the value of the specified keywords versus timestep.
        :param keyword: name of keywords in LAMMPS log file, e.g. TotEng, Temp, Volume
        :return: None
        """

        new_df = self.df[keyword].copy()
        if len(keyword) == 1:
            ax = new_df.plot()
        elif len(keyword) == 2:
            ax = new_df.plot(secondary_y=keyword[1:], legend=True)
        else:
            ax = new_df.plot()
        ax.set_xlabel("Timestep")

        plt.show()

    def save_txt(self, *keyword):
        new_df = self.df[list(keyword)].copy()

    def estimate_time(self, timestep=0):
        if not timestep:
            return

        line = np.polyfit(self.df.index, self.df['CPU'], 1)
        est_time = line[0] * timestep + line[1]
        print("Estimated time: about %s seconds from the start" % int(est_time))
        print("Estimated end time: about %s seconds later" % int(est_time - list(self.df['CPU'])[-1]))
