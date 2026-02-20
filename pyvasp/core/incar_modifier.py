class IncarModifier:

    def __init__(self, job, cluster=None):
        self.job = job
        self.cluster = cluster
        self.change = {}
        self.remove = []

    def apply_job(self):
        rule = self.job.get_rules()
        self.change.update(rule.get("change", {}))
        self.change.update(rule.get("active", {}))
        self.remove.extend(rule.get("out", []))

    def apply_cluster(self):
        if not self.cluster:
            return
        rule = self.cluster.get_incar_rules()
        self.change.update(rule.get("change", {}))

    def apply_cli(self, cli_dict):
        if cli_dict:
            self.change.update(cli_dict)

    def build(self):
        self.apply_job()
        self.apply_cluster()
        return self.change, self.remove
