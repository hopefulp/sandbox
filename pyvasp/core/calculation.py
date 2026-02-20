from .incar_modifier import IncarModifier
from vaspflow.io.incar_io import modify_incar

class VaspCalculation:

    def __init__(self, job, cluster, queue=None):
        self.job = job
        self.cluster = cluster
        self.queue = queue
        self.incar_change = {}
        self.incar_remove = []

    def configure_incar(self, cli_override=None):
        modifier = IncarModifier(self.job, self.cluster)
        modifier.apply_cli(cli_override)
        self.incar_change, self.incar_remove = modifier.build()

    def write_incar(self, path="INCAR"):
        modify_incar(path, self.incar_change, self.incar_remove)

    def submit(self):
        cmd = self.cluster.build_submit_command(
            jobname=self.job.name,
            queue=self.queue
        )
        print(cmd)
