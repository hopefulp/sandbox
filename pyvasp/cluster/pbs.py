from .base import Cluster

class PbsCluster(Cluster):

    def build_submit_command(self, jobname, queue):
        return f"qsub -N {jobname} pbs_vasp.sh"
