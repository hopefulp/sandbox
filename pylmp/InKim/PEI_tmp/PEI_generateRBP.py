#!/usr/bin/python3

"""
Generate N random hyperbranched polymers of n monomers by random addition method.
"""

import sys
import getopt
import glob
import time
import pickle
import matplotlib
# this Agg mode enables me to use matplotlib without X display: cool!
# https://crazythinking.wordpress.com/2011/05/01/using-matplotlib-pylab-without-a-display/
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from monomer import *


version = "20170220"


def main(ter, lin, den, n_sample, ff, n_monomer, prefix, ratio=0, compile=False, silent=False):

    # workspace init
    curr_dir = os.path.abspath("./")
    out_dir = curr_dir + "/" + prefix + "/"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # log init
    log_filename = "generateRBP.log"
    if not silent:
        print("log file will be: " + log_filename)
    log_file = open(log_filename, 'a', buffering=1)
    output = "Date\tFilename\tN_monomer\tDB\tWI\tDirectory\tElapsed(s)\n"
    output += "Terminal unit: %s\n" % os.path.abspath(ter)
    output += "Linear unit: %s\n" % os.path.abspath(lin)
    output += "Dendron unit: %s\n" % os.path.abspath(den)

    log_file.write(output)

    # avoid file overwrite
    counter = 0
    filename_tester = out_dir + prefix + "_"
    pool = glob.glob(filename_tester + "*.bgf")
    if len(pool) > 0:
        j = [i.split(filename_tester)[-1] for i in pool]
        j = [i.strip(".bgf") for i in j]
        j = [int(i) for i in j]
        j.sort()    # j[-1] has the maximum
        counter = j[-1]

    # store WI index
    wi_pkl_file = prefix + ".WI.pickle"
    wi_file = prefix + ".WI.dat"
    if os.path.exists(wi_pkl_file):
        with open(wi_pkl_file, 'rb') as f:
            wis = pickle.load(f)
    else:
        wis = []

    # generation
    i = 0
    trial = 0
    while (i < n_sample):
        trial += 1
        if i != 0 and trial > (i * 10):
            print("STOP: Exceeding maximum number of trials.")
            break;
        filename = out_dir + prefix + "_" + str(i + counter)    # first guess
        jobname = prefix + "_" + str(i + counter)
        if not silent:
            print("\n** Generating %s (trial: %s)" % (jobname, trial))

        # random polymer init
        g = RandomPolymer()
        g.Dendron = den
        g.Linear = lin
        g.Terminal = ter
        g.ff = ff
        g.num_target_node = n_monomer - 1

        # random polymer generation 
        g.add_focal_point()
        g.build_random()

        # Degree of Branching and Wiener Index
        db = g.calculate_db()
        if ratio and not (ratio * 0.9 <= db <= ratio * 1.1):
            print("DISCARD: DB not within of 10% of the target DB.")
            continue

        wi = g.calculate_wiener_index()
        if wi in wis:
            print("DISCARD: Same WI structure found.")
            continue
        else:
            wis.append(wi)

        # draw connection graph
        # edge label
        lbl = dict()
        for j in g.edges():
            # lbl[i] = str(i[0]) + "_" + str(g[i[0]][i[1]]['branch'])
            lbl[j] = str(g[j[0]][j[1]]['branch'])

        # decoration
        pos = nx.spring_layout(g, k=0.5)
        nx.draw_networkx_edges(g, pos)
        nx.draw_networkx_labels(g, pos)
        nx.draw_networkx_edge_labels(g, pos, edge_labels=lbl, label_pos=0.7, font_size=9)
        g.set_nodes_degree_of_branching()
        t_nodes = [i for (i, j) in g.nodes(data=True) if g[i]['n_branch'] == 0]
        l_nodes = [i for (i, j) in g.nodes(data=True) if g[i]['n_branch'] == 1]
        d_nodes = [i for (i, j) in g.nodes(data=True) if g[i]['n_branch'] == 2]
        nx.draw_networkx_nodes(g, pos, nodelist=t_nodes, node_color='w')
        nx.draw_networkx_nodes(g, pos, nodelist=l_nodes, node_color='y')
        nx.draw_networkx_nodes(g, pos, nodelist=d_nodes, node_color='g')
        nx.draw_networkx_nodes(g, pos, nodelist=[0], node_color='r')
        g.remove_nodes_degree_of_branching()

        # connection graph save as PNG file
        plt.savefig(filename + '.png')
        if not silent:
            print("\tFigure saved as PNG file: " + filename + ".png")

        # save graph
        nx.write_gpickle(g, filename + ".gpickle")
        if not silent:
            print("\tGraph saved as gpickle file: " + filename + ".gpickle")

        # compile structure and save
        if compile:
            t1 = time.time()
            if not silent:
                print("\tCompiling..")
            g.compile(filename + ".bgf")
            if not silent:
                print("\tModel saved as BGF file: " + filename + ".bgf")
            t2 = time.time()
            elapsed_time = t2 - t1

            # calculate Mw
            mw = bgftools.getMass(g.bgfmodel, ff_file=g.ff)
        else:
            elapsed_time = 0.0
            mw = 0

        # write log
        output = str(time.asctime(time.gmtime())) + "\t" + str(jobname) + "\t" + str(n_monomer) + "\t" + "{0:8.3f}".format(mw) + "\t" + "{0:8.3f}".format(db) + "\t" + str(int(wi)) + "\t" + str(filename) + "\t" + str(elapsed_time) + "\n"
        log_file.write(output)

        if not silent:
            print(">> Generated a polymer with " + str(len(g.nodes())) + " monomers..")

        with open(wi_file, mode='a') as f:
            f.write("%s %s\n" % (filename, wi))

        i += 1

    # record WI
    with open(wi_pkl_file, 'wb') as f:
        pickle.dump(wis, f)

    # end of main()


if __name__ == "__main__":
    option = ""
    args = ""
    is_compile = False
    usage = "Version " + str(version) + """

Usage: PEI_generateRBP.py OPTIONS

Generates random hyperbranched polymers by slow addition method with Degree of Branching "unspecified."
If you want to specify DB in terms of amine ratio, use PEI_generate.py instead.

    Required:
        -t terminal unit in BGF
        -l linear unit in BGF
        -d dendron unit in BGF
        -N number of polymer generation (Default: 1)
        -n number of monomers in a polymer
        -D output_directory (Default: .)
        -f ff_file (Default: DREIDING 2)
        -s prefix (Default: randomPolymer)
        -c compile? (Default: yes)

    * Monomer units must have chain H for head atom and chain T for branch atom
    * Polymer generation list will be recorded to the file {prefix}.log

Please report any bugs to in.kim@kaist.ac.kr

    """

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    # Defaults
    silent = False
    T = L = D = ff = prefix = ""
    N = 0   # number of generating polymers
    n = 0   # number of monomers
    r = 0.0 # degree of branching: DB = 2D/(T+L+D)

    options, args = getopt.getopt(sys.argv[1:], 'ht:l:d:N:f:n:s:cr:',
                                  ['help', 'terminal=', 'linear=', 'dendron=', 'number=', 'ff=', 'n=', 'prefix=', 'compile', 'ratio='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0)
        elif option in ('-t', '--terminal'):
            T = os.path.abspath(value)
        elif option in ('-l', '--linear'):
            L = os.path.abspath(value)
        elif option in ('-d', '--dendron'):
            D = os.path.abspath(value)
        elif option in ('-N', '--number'):
            N = int(value)
        elif option in ('-s', '--prefix'):
            prefix = str(value)
        elif option in ('-f', '--ff'):
            ff = str(value)
        elif option in ('-n', '--nmonomer'):
            n = int(value)
        elif option in ('-c', '--compile'):
            is_compile = True
        elif option in ('-r', '--ratio'):
            r = float(value)
        else:
            print(usage)
            sys.exit(0)

    # defaults
    if ff == "":
        ff = "/home/noische/ff/DREIDING2.21.ff"
    if N == 0:
        N = 1
    if n == 0:
        nu.die("You must specify number of monomers in a hb polymer.")
    if prefix == "":
        prefix = "randomPolymer"

    # main run
    main(T, L, D, N, ff, n, prefix, ratio=r, compile=is_compile, silent=silent)
