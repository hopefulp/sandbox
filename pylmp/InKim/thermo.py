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

# global variables
extensive_properties = ['natom', 'nmolecules', 'temperature', 'volume', 'pressure']
extensive_properties_title = ['natom', 'nmol', 'temp', 'vol', 'press']
extensive_properties_format = ['%12d', '%12.3f', '%12.3f', '%12.3f', '%12.3f']
twopt_properties = ['dof', 'fluidicity', 'Diffus']
twopt_properties_title = ['dof', 'f', 'Df']
twopt_properties_format = ['%14.4f', '%14.6e', '%14.6e']
basic_format = "%8.3f"

class thermofile():
    def __init__(self, thermo_file):
        # class init
        self.properties = []
        self.thermo = nu.hash() # e.g) thermo[group_no]['Trns']['Sq']
        self.thermo_per_mol = nu.hash()
        self.timestep = -1

        # local init
        pat_keyw = re.compile(r'(\S+)\s*\(\S+\)')
        pat_cat = re.compile(r'(\S+)\[\w(\d+)\]')
        n_group = 0
        _header = ""

        # read
        #print("Reading thermo file %s" % thermo_file)

        try:
            self.timestep = int(re.search(re.compile(r'\.(\d+)\.2pt\.mol\.grps\.thermo'), thermo_file).group(1))
        except:
            nu.warn('No timestep found in the filename %s' % thermo_file)
            return

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


#def main(directory, nsamples, out_file, print_stdev=False, report=False, silent=True):
class thermo():
    def __init__(self, directory):
        self.data = nu.hash()
        self.categories = []
        self.properties = []

        # init
        self.data = nu.hash()   # data[timestep][group][category][property]
        #REMARK: categories = ["Tgas", "Tsol", "Trns", "Rgas", "Rsol", "Rot", "Ivib", "Tot"]
        #REMARK: properties = ["temperature", "pressure", "volume", "ZPE", "Emd", "Eq", "Sq", "Aq", "Cvq", "S(0)", "Diffus", "fluidicity"]
        
        self.path = glob.glob(directory + "/*thermo")

        for f in self.path:
            _ = thermofile(f)
            self.data[_.timestep] = _.thermo_per_mol

        self.timesteps = sorted(self.data.keys())
        self.groups = self.data[self.timesteps[0]].keys()
        self.categories = self.data[self.timesteps[0]][self.groups[0]].keys()
        self.properties = self.data[self.timesteps[0]][self.groups[0]][self.categories[0]].keys()


    def evolution(self, group, category, property):
        """
        print time evolution of values
        """
        output = "t\t%s_%s_%s\n" % (group, category, property)
        for t in self.timesteps:
            output += "{:<12d}".format(t) + basic_format % self.data[t][group][category][property] + "\n"

        print(output)


    def do_average(self, nsamples=0):
        """
        calculate average and stdev over timestep
        """
        self.sum = nu.hash(); self.sum_sq = nu.hash(); self.average = nu.hash(); self.stdev = nu.hash();

        if nsamples:
            self.samples = self.timesteps[-nsamples:]
        else:
            self.samples = self.timesteps

        print(self.samples)

        n = float(len(self.samples))

        for t in self.data:
            if t in self.samples:
                for g in self.data[t]:
                    for c in self.data[t][g]:
                        for p in self.data[t][g][c]:
                            if p in self.sum[g][c]:
                                self.sum[g][c][p] += self.data[t][g][c][p]
                                self.sum_sq[g][c][p] += self.data[t][g][c][p]**2
                            else:
                                self.sum[g][c][p] = 0.0
                                self.sum_sq[g][c][p] = 0.0
                                self.average[g][c][p] = 0.0
                                self.stdev[g][c][p] = 0.0
                                self.sum[g][c][p] += self.data[t][g][c][p]
                                self.sum_sq[g][c][p] += self.data[t][g][c][p]**2

        for g in self.sum:
            for c in self.sum[g]:
                for p in self.sum[g][c]:
                    self.average[g][c][p] = self.sum[g][c][p] / n
                    self.stdev[g][c][p] = (self.sum_sq[g][c][p] - self.sum[g][c][p]**2 / n ) / n


    def print_data(self, hash):
        '''
        print sum or averaged data
        '''
        '''
        # print
        g = hash.keys()
        for i in g:
            print("Group %d" % i)
            print(" " * 8 + "".join("%12s" % i for i in extensive_properties_title) + "".join("%12s" % i for i in properties) + "".join('%14s' % i for i in twopt_properties_title))
            for c in categories:
                print("%8s" % c + "".join(extensive_properties_format[index] % hash[i][c][p] for index, p in enumerate(extensive_properties)) + "".join("%12.4f" % hash[i][c][p] for p in properties) + "".join(twopt_properties_format[index] % hash[i][c][p] for index, p in enumerate(twopt_properties)))
            print("")

        '''
