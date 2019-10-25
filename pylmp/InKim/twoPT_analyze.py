#!/home/noische/python

import sys
import glob
import pickle as pkl
import re
import copy
import getopt
import numpy as np
import nutils as nu

version = '160502'

class thermo():

    def __init__(self, thermo_file):
        # class init
        self.properties = []
        self.thermo = nu.hash() # e.g) thermo[group_no]['Trns']['Sq']
        self.thermo_per_mol = nu.hash()
        self.timestep = -1

        # local init
        pat_keyw = re.compile(r'(\S+)\s*\(\S+\)')
        pat_cat = re.compile(r'(\S+)\[\w(\d+)\]')
        categorys = []  # not a typo, it's a personal onvention!
        groups = []
        n_group = 0
        _header = ""

        # read
        #print("Reading thermo file %s" % thermo_file)

        try:
            self.timestep = int(re.search(re.compile(r'\.(\d+)\.2pt\.mol\.grps\.thermo'), thermo_file).group(1))
        except:
            nu.warn('No timestep found in the filename %s' % thermo_file)
            return 0

        with open(thermo_file, 'r') as f:
            for line in f:
                line = line.replace('\n', '')

                if "*" in line:
                    continue
                if "Calculation" in line:
                    continue

                if line:
                    parse = " ".join(line.split())
                    parse = parse.split()

                    # header
                    if 'property' in parse:
                        self._header = parse
                    else:
                        property = parse[0].replace('_', '')
                        property = property.replace('/', '')
                        if re.search(pat_keyw, property):
                            property = re.search(pat_keyw, property).group(1)

                        self.properties.append(property)

                        for index, i in enumerate(parse[1:]):
                            _ = self._header[index+1]  # e.g) Trns[G001]
                            cat = re.search(pat_cat, _).group(1)    # Trns
                            grp = int(re.search(pat_cat, _).group(2))   # 1

                            self.thermo[grp][cat][property] = float(i)


        # per molecule values
        self.thermo_per_mol = copy.deepcopy(self.thermo)
        skip_properties = ['nmolecules', 'natom', 'temperature', 'pressure', 'volume', 'Diffus', 'fluidicity', 'dof']
        for g in self.thermo.keys():
            for c in self.thermo[g].keys():
                for p in self.thermo[g][c].keys():
                    if not p in skip_properties:
                        self.thermo_per_mol[g][c][p] /= self.thermo_per_mol[g][c]['nmolecules']


def main(directory, nsamples, out_file, print_stdev=False, report=False, silent=True):
    # init
    d = nu.hash()   # d[timestep][group][category][property]
    #categories = ["Tgas", "Tsol", "Trns", "Rgas", "Rsol", "Rot", "Ivib", "Tot"]
    #properties = ["temperature", "pressure", "volume", "ZPE", "Emd", "Eq", "Sq", "Aq", "Cvq", "S(0)", "Diffus", "fluidicity"]
    
    # glob *thermo
    thermos = glob.glob(directory + "/*thermo")

    # read 2pt data file
    for f in thermos:
        _ = thermo(f)
        d[_.timestep] = _.thermo_per_mol

    # determine #timestep to average
    timesteps = list(d.keys())
    timesteps.sort()
    if nsamples:
        t_sample = timesteps[-nsamples:]
    else:
        t_sample = timesteps

    print(t_sample)
    print("")

    # create report
    t0 = t_sample[0]
    groups = list(d[t0].keys())
    categories = list(d[t0][groups[0]].keys())
    properties = list(d[t0][groups[0]][categories[0]].keys())
    extensive_properties = ['natom', 'nmolecules', 'temperature', 'volume', 'pressure']
    extensive_properties_title = ['natom', 'nmol', 'temp', 'vol', 'press']
    extensive_properties_format = ['%12d', '%12.3f', '%12.3f', '%12.3f', '%12.3f']
    twopt_properties = ['dof', 'fluidicity', 'Diffus']
    twopt_properties_title = ['dof', 'f', 'Df']
    twopt_properties_format = ['%14.4f', '%14.6e', '%14.6e']

    for i in extensive_properties:
        properties.remove(i)
    for i in twopt_properties:
        properties.remove(i)

    if report:
        important_properties = ["Sq", "Aq", "Eq"]
        important_categories = ["Trns", "Rot", "Ivib", "Tot"]
        for g in d[t0]:
            log_file = "2pt.report.group" + "%s" % g + ".log"
            with open(log_file, 'w') as f:
                output = "t\t" + "\t".join("%s_%s" % (c, p) for c in important_categories for p in important_properties) + "\n"
                for t in sorted(d.keys()):
                    output += "{:<12d}".format(t) + " ".join("%8.3f" % d[t][g][c][p] for c in important_categories for p in important_properties) + "\n"
                f.write(output)
            print("* A file %s is generated for group %s." % (log_file, g))
    print("")

    # calculate average and stdev over timestep
    sum = nu.hash(); sum_sq = nu.hash(); average = nu.hash(); stdev = nu.hash();

    n = float(len(t_sample))
    for t in d:
        if t in t_sample:
            for g in d[t]:
                for c in d[t][g]:
                    for p in d[t][g][c]:
                        if p in sum[g][c]:
                            sum[g][c][p] += d[t][g][c][p]
                            sum_sq[g][c][p] += d[t][g][c][p]**2
                        else:
                            sum[g][c][p] = 0.0
                            sum_sq[g][c][p] = 0.0
                            average[g][c][p] = 0.0
                            stdev[g][c][p] = 0.0
                            sum[g][c][p] += d[t][g][c][p]
                            sum_sq[g][c][p] += d[t][g][c][p]**2

    for g in sum:
        for c in sum[g]:
            for p in sum[g][c]:
                average[g][c][p] = sum[g][c][p] / n
                stdev[g][c][p] = (sum_sq[g][c][p] - sum[g][c][p]**2 / n ) / n

    # print
    print("** Data averaged over %d snapshots **" % int(n))
    if print_stdev:
        print("** Stdev **")
        g = stdev.keys()
        for i in g:
            print("Group %d" % i)
            print(" " * 8 + "".join("%12s" % i for i in extensive_properties_title) + "".join("%12s" % i for i in properties) + "".join('%14s' % i for i in twopt_properties_title))
            for c in categories:
                #print("" + "%8s" % c + "  " + "".join("%12.4f" % stdev[i][c][p] for p in properties))
                print("%8s" % c + "".join(extensive_properties_format[index] % stdev[i][c][p] for index, p in enumerate(extensive_properties)) + "".join("%12.4f" % stdev[i][c][p] for p in properties) + "".join(twopt_properties_format[index] % stdev[i][c][p] for index, p in enumerate(twopt_properties)))
            print("")

    else:
        print("** Mean **")
        g = average.keys()
        for i in g:
            print("Group %d" % i)
            print(" " * 8 + "".join("%12s" % i for i in extensive_properties_title) + "".join("%12s" % i for i in properties) + "".join('%14s' % i for i in twopt_properties_title))
            for c in categories:
                print("%8s" % c + "".join(extensive_properties_format[index] % average[i][c][p] for index, p in enumerate(extensive_properties)) + "".join("%12.4f" % average[i][c][p] for p in properties) + "".join(twopt_properties_format[index] % average[i][c][p] for index, p in enumerate(twopt_properties)))
            print("")

    print("Done")
    

if __name__ == "__main__":

    # Initialize
    n = 0;
    directory = "results"
    out_file = "2pt_analyze.profile"
    use_stdev = False
    print_report = False
    usage = """
    """

    print("\n" + sys.argv[0] + " version " + str(version) + "\n")

    if len(sys.argv) < 2:
        print(usage);
        sys.exit(1)

    options, args = getopt.getopt(sys.argv[1:], 'hd:n:o:sr', ['help','directory=','nsamples=','out=', 'stdev', 'report'])
    try:
        for option, value in options:
            if option in ('-h', '--help'):
                print(usage)
                sys.exit(0);
            elif option in ('-d', '--directory'):
                directory = value
            elif option in ('-n', '--nsamples'):
                n = int(value)
            elif option in ('-o', '--out'):
                out_file = value
            elif option in ('-s'):
                use_stdev = True
            elif option in ('-r'):
                print_report = True
            elif not option:
                print(usage)
                sys.exit(0)
    except ValueError:
        nu.die("script cannot continue. Check input parameters.")

    main(directory, n, out_file, print_stdev=use_stdev, report=print_report, silent=True)
