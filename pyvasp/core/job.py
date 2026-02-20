from .incar_rules import JOB_RULES

class Job:

    def __init__(self, name):
        self.name = name

    def get_rules(self):
        return JOB_RULES.get(self.name, {})
