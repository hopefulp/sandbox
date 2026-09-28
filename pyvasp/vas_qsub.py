#!/home/joonho/anaconda3/bin/python

from dataclasses import dataclass
from typing import Optional

from libcluster import detect_cluster
from server_env import nXn
from common import yes_or_no
from parsing import startnum

import os
import sys


# ==========================================================
# Data model
# ==========================================================

@dataclass
class QueueConfig:
    partition: int        # e.g. 1~6
    nnode: int
    nproc: Optional[int] = None


# ==========================================================
# Public API
# ==========================================================

def qsub_command(
    ndir,
    queue: Optional[QueueConfig] = None,
    option: Optional[str] = None,
    vasp_exe: Optional[str] = None,
    lkisti: Optional[str] = None,
    Lrun: Optional[bool] = None,
    cluster: Optional[str] = None,
):
    """
    Build and optionally execute job submission command.
    """

    cluster = cluster or detect_cluster()

    if cluster == "kisti":
        cmd = _build_kisti_command(ndir, queue, option, vasp_exe, lkisti)

    elif cluster == "pt":
        if not queue:
            raise ValueError("QueueConfig required for pt cluster")

        cmd = _build_slurm_command(ndir, queue, option, vasp_exe)

    else:
        raise RuntimeError("No submission rule for cluster: {}".format(cluster))

    print(cmd)

    if Lrun or yes_or_no("Will you run in qsub?"):
        os.system(cmd)

    return cmd


# ==========================================================
# KISTI (PBS)
# ==========================================================

def _build_kisti_command(ndir, queue, option, vasp_exe, lkisti):

    str_vasp = _kisti_vasp_flag(vasp_exe)

    nnode = 20
    np = 40

    # A supplied queue overrides the resource line embedded in the PBS script.
    # NEB supplies intermediate_images * nodes_per_image.
    if queue:
        nnode = queue.nnode
        np = queue.nproc or np
        vasp_flag = "{} ".format(str_vasp) if str_vasp else ""
        return (
            "qsub -N {} {}"
            "-l select={}:ncpus=40:mpiprocs={}:ompthreads=1 "
            "$SB/pypbs/pbs_vasp_kisti_skl.sh"
        ).format(ndir, vasp_flag, nnode, np)

    if option == "mem":
        hproc = np // 2
        return (
            "qsub -N {} "
            "-l select={}:ncpus={}:mpiprocs={}:ompthreads=1 "
            "$SB/pypbs/pbs_vasp_kisti_skl.sh"
        ).format(ndir, nnode, np, hproc)

    if option and "long" in option:
        hour = _extract_hour(option, default=96)
        return (
            "qsub -N {} -q long "
            "-l walltime={}:00:00 "
            "$SB/pypbs/pbs_vasp_kisti_skl.sh"
        ).format(ndir, hour)

    if option == "opt":
        return "qsub -N {} $SB/pypbs/pbs_vasp_kisti_sklopt.sh".format(ndir)

    if lkisti == "kp":
        return (
            "qsub -N {} {} "
            "-l walltime=1:00:00 "
            "$SB/pypbs/pbs_vasp_kisti_skl.sh"
        ).format(ndir, str_vasp)

    return "qsub -N {} {} $SB/pypbs/pbs_vasp_kisti_skl.sh".format(ndir, str_vasp)


def _kisti_vasp_flag(vasp_exe):
    if not vasp_exe:
        return ""

    if "g" in vasp_exe:
        return "-v exe=gam"
    if "xy" in vasp_exe:
        return "-v exe=xyrelax"
    if "ncl" in vasp_exe:
        return "-v exe=ncl"

    return ""


# ==========================================================
# PT (SLURM)
# ==========================================================

def _build_slurm_command(ndir, queue, option, vasp_exe):

    X = queue.partition
    nnode = queue.nnode

    # total MPI processes
    if queue.nproc:
        nproc = queue.nproc
    else:
        nproc = nnode * nXn[X]

    str_vasp = _slurm_vasp_flag(vasp_exe)

    if option == "mem":
        hproc = nXn[X] // 2
        return (
            "sbatch -J {} -p X{} "
            "-N {} -c {} "
            "--export=hmem=1 "
            "/home/joonho/sandbox/pypbs/slurm_sbatch.sh"
        ).format(ndir, X, nnode, hproc)

    if option == "opt":
        return (
            "sbatch -J {} -p X{} "
            "-N {} -n {} "
            "{} "
            "/home/joonho/sandbox/pypbs/slurm_sbatch_vaspopt.sh"
        ).format(ndir, X, nnode, nproc, str_vasp)

    if option == "sim":
        return (
            "sbatch -J {} -p X{} "
            "-N {} -n {} "
            "{} "
            "/home/joonho/sandbox/pypbs/slurm_sbatch_sim.sh"
        ).format(ndir, X, nnode, nproc, str_vasp)

    return (
        "sbatch -J {} -p X{} "
        "-N {} -n {} "
        "{} "
        "/home/joonho/sandbox/pypbs/slurm_sbatch.sh"
    ).format(ndir, X, nnode, nproc, str_vasp)


def _slurm_vasp_flag(vasp_exe):
    if not vasp_exe:
        return ""

    if "g" in vasp_exe:
        return "--export=exe=gam"
    if "xy" in vasp_exe:
        return "--export=exe=xyrelax"
    if "ncl" in vasp_exe:
        return "--export=exe=ncl"

    return ""


# ==========================================================
# Utilities
# ==========================================================

def _extract_hour(option, default=96):
    if option.isalpha():
        return default

    idx = startnum(option)
    return option[idx:] if idx is not None else default
