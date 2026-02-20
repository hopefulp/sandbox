#!/home/joonho/anaconda3/bin/python

from dataclasses import dataclass
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
    partition: int   # e.g. 1~6
    nnode: int
    nproc: int | None = None


# ==========================================================
# Public API
# ==========================================================

def qsub_command(
    ndir,
    queue: QueueConfig | None = None,
    option: str | None = None,
    vasp_exe: str | None = None,
    lkisti: str | None = None,
    Lrun: bool | None = None,
    cluster: str | None = None,
):
    """
    Build and optionally execute job submission command.
    """

    cluster = cluster or detect_cluster()

    if cluster == "kisti":
        cmd = _build_kisti_command(ndir, option, vasp_exe, lkisti)

    elif cluster == "pt":
        if not queue:
            raise ValueError("QueueConfig required for pt cluster")

        cmd = _build_slurm_command(ndir, queue, option, vasp_exe)

    else:
        raise RuntimeError(f"No submission rule for cluster: {cluster}")

    print(cmd)

    if Lrun or yes_or_no("Will you run in qsub?"):
        os.system(cmd)

    return cmd


# ==========================================================
# KISTI (PBS)
# ==========================================================

def _build_kisti_command(ndir, option, vasp_exe, lkisti):

    str_vasp = _kisti_vasp_flag(vasp_exe)

    nnode = 20
    np = 40

    if option == "mem":
        hproc = np // 2
        return (
            f"qsub -N {ndir} "
            f"-l select={nnode}:ncpus={np}:mpiprocs={hproc}:ompthreads=1 "
            f"$SB/pypbs/pbs_vasp_kisti_skl.sh"
        )

    if option and "long" in option:
        hour = _extract_hour(option, default=96)
        return (
            f"qsub -N {ndir} -q long "
            f"-l walltime={hour}:00:00 "
            f"$SB/pypbs/pbs_vasp_kisti_skl.sh"
        )

    if option == "opt":
        return f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_sklopt.sh"

    if lkisti == "kp":
        return (
            f"qsub -N {ndir} {str_vasp} "
            f"-l walltime=1:00:00 "
            f"$SB/pypbs/pbs_vasp_kisti_skl.sh"
        )

    return f"qsub -N {ndir} {str_vasp} $SB/pypbs/pbs_vasp_kisti_skl.sh"


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

def _build_slurm_command(ndir, queue: QueueConfig, option, vasp_exe):

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
            f"sbatch -J {ndir} -p X{X} "
            f"-N {nnode} -c {hproc} "
            f"--export=hmem=1 "
            f"/home/joonho/sandbox/pypbs/slurm_sbatch.sh"
        )

    if option == "opt":
        return (
            f"sbatch -J {ndir} -p X{X} "
            f"-N {nnode} -n {nproc} "
            f"{str_vasp} "
            f"/home/joonho/sandbox/pypbs/slurm_sbatch_vaspopt.sh"
        )

    if option == "sim":
        return (
            f"sbatch -J {ndir} -p X{X} "
            f"-N {nnode} -n {nproc} "
            f"{str_vasp} "
            f"/home/joonho/sandbox/pypbs/slurm_sbatch_sim.sh"
        )

    return (
        f"sbatch -J {ndir} -p X{X} "
        f"-N {nnode} -n {nproc} "
        f"{str_vasp} "
        f"/home/joonho/sandbox/pypbs/slurm_sbatch.sh"
    )


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
