from .slurm import SlurmCluster
from .pbs import PbsCluster
from vaspflow.infra.detect import detect_cluster

def create_cluster():

    name = detect_cluster()

    if name == "pt":
        return SlurmCluster()
    elif name == "kisti":
        return PbsCluster()
    else:
        raise RuntimeError("Unknown cluster")
