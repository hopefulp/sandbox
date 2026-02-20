from .base import Cluster

class SlurmCluster(Cluster):

    def build_submit_command(self, jobname, queue):

        X = queue.partition
        nnode = queue.nnode
        nproc = queue.nproc or 32

        return (
            f"sbatch -J {jobname} "
            f"-p X{X} "
            f"-N {nnode} "
            f"-n {nproc} "
            "slurm_sbatch.sh"
        )
