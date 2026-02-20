from abc import ABC, abstractmethod

class Cluster(ABC):

    @abstractmethod
    def build_submit_command(self, jobname, queue):
        pass

    def get_incar_rules(self):
        return {}
