#!/home/joonho/anaconda3/envs/py310/bin/python

import argparse
from vaspflow.core.job import Job
from vaspflow.core.calculation import VaspCalculation
from vaspflow.cluster.factory import create_cluster
from vaspflow.infra.queue import QueueConfig

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--job", required=True)
    parser.add_argument("-p", "--partition", type=int)
    parser.add_argument("-N", "--nnode", type=int)
    parser.add_argument("-np", "--nproc", type=int)
    parser.add_argument("--incar", nargs="*", help="override INCAR KEY=VAL")

    args = parser.parse_args()

    job = Job(args.job)
    cluster = create_cluster()

    queue = None
    if args.partition and args.nnode:
        queue = QueueConfig(args.partition, args.nnode, args.nproc)

    cli_dict = {}
    if args.incar:
        for kv in args.incar:
            k, v = kv.split("=")
            cli_dict[k] = v

    calc = VaspCalculation(job, cluster, queue)
    calc.configure_incar(cli_override=cli_dict)
    calc.write_incar()
    calc.submit()


if __name__ == "__main__":
    main()
