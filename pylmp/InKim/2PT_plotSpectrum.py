#!/home/noische/Enthought/Canopy_64bit/User/bin/python

"""
2PT_plotSpectrum.py: Read power spectrum of 2PT results (*.pwr) and interactively plot.

** Rules **
1. timesteps are in the filename.

"""
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

import sys
import os
import math
import pickle
import glob
import re
import nutils as nu


curr_t = 0; curr_group = 0; curr_pwr_type = '';

def parse_pwr_files():

    pwr_files = nu.hash();
    data = nu.hash();
    pwr_types = ["PWR_hs", "PWR_st", "PWR_fr", "PWR_sr", "PWRcmt", "PWRcmo", "PWRimv", "PWRcmr", "PWRtot"]
    pwr_type_pat = re.compile(r'([A-z]+)\[\w(\d+)\]')

    # Read timesteps
    l_timestep = []
    pwr_filenames = glob.glob("*2pt.mol.grps.pwr")
    timestep_pat = re.compile(r"\.([0-9]*)\.")

    for fname in pwr_filenames:
        timestep = int(re.findall(timestep_pat, fname)[0])
        pwr_files[timestep]['filename'] = fname
    #pwr_files[0]['filename'] = 'tester2.pwr'    #### test

    print("%d 2PT power spectrum files found." % len(pwr_files))

    # Read power spectrum data in a file
    for i, timestep in enumerate(pwr_files):
        #if i == 5:
        #    break   # tester

        _ = pwr_files[timestep]['filename']
        f = open(pwr_files[timestep]['filename'])
        sys.stdout.write("\r  ..reading a file %s .. %d/%d" % (_, i, len(pwr_files)))
        sys.stdout.flush()
        cols, col_index = nu.get_columns(f, delim=" ", header=3, type="float")

        # Generate a dictionary
        """
        data[timestep][pwr_type][group] = [values.....]
        data[timestep]['freq(cm-1)'] = [0.0, 1.1119, 2.2238, 3.3356, 4.4475, 5.5594, 6.6713, 7.7832, 8.895, 10.0069, ...]
        pwr_type: PWR_hs, PWR_st, PWR_fr, PWR_sr, PWRcmt, PWRcmo, PWRimv, PWRcmr, PWRtot
        """
        for col in cols:
            try:
                pwr_type, group = re.findall(pwr_type_pat, col)[0]
            except IndexError:
                data[timestep]["freq(cm-1)"] = cols[col]
            else:
                data[timestep][pwr_type][int(group)] = cols[col]

        i += 1
    
    # Returns the dictionary
    return data


def main():

    ### initialize
    fig = plt.figure()
    ax = fig.add_subplot(111)
    axcolor = 'lightgoldenrodyellow'

    ### Read data
    if len(sys.argv) == 1:
        print("No spectrum pickle file given. The script will extract and write data to pwr_spectrum.pickle")
        result = parse_pwr_files()
        f = open("pwr_spectrum.pickle", 'w')
        pickle.dump(result, f)
        f.close()
        print("Generated: pwr_spectrum.pickle")
    else:
        try:
            datafile = sys.argv[1]
            print("Reading pickle file %s" % datafile)
            f = open(datafile)
            result = pickle.load(f)
            f.close()
        except:
            nu.die("Please specify a data or run the script without any parameters.")

    ### REMARK: there's a problem that n and r cannot be loaded in the eventhandler.
    ###         so add a global value and set when the r changes.
    # power_pwr_types
    def get_pwr_type():
        global curr_pwr_type
        return curr_pwr_type
    def set_pwr_type(val):
        global curr_pwr_type
        curr_pwr_type = val

    # timesteps
    def get_t():
        global curr_t
        return curr_t
    def set_t(val):
        global curr_t
        curr_t = int(val)

    # groups 
    def get_group():
        global curr_group
        return curr_group
    def set_group(val):
        global curr_group
        curr_group = int(val)

    ### get the data from the result dictionary
    # result: [timestep][pwr_type][group]
    timesteps = result.keys()
    timesteps.sort()
    dict_t = {}         # dict_t[i] = timestep: bookkeeping
    for index, value in enumerate(timesteps):
        dict_t[index] = value 
    #print timesteps
    #print dict_t
    t = timesteps[0]	# the first timestep
    set_t(t)
    n_t = len(timesteps) - 1	# number of timesteps

    pwr_types = result[t].keys()
    for i in pwr_types:
        if "freq" in i:
            pwr_types.remove(i)

    pwr_type = pwr_types[0]
    set_pwr_type(pwr_type)

    groups = result[t][pwr_type].keys()
    group = groups[0]
    set_group(group)

    spectrum = result[t][pwr_type][group]

    x = result[t]['freq(cm-1)']

    '''
    plt.autoscale(True)
    l, = plt.plot(x, spectrum)
    ax = plt.gca()
    #ax.autoscale_view(True, True, True)
    ax.set_xlim([0, 1500.0])
    ax.set_ylim([0, 30.0])


    ### timestep slider
    ax_timestep = plt.axes([0.25, 0.03, 0.65, 0.03], axisbg=axcolor)
    #s_timestep = Slider(ax_timestep, 'Timestep', 0, n_t, valinit=t, valfmt='%8d')
    s_timestep = Slider(ax_timestep, 'Timestep', 0, n_t, valinit=dict_t[0], valfmt='%4.0f', dragging=False)

    def t_update(val):
        g = get_group();
        pwr_type = get_pwr_type();

        t = dict_t[round(val)]
        s_timestep.valtext.set_text(str(t))
        set_t(t)

        print t, pwr_type, g

        l.set_ydata(result[t][pwr_type][g])
        plt.draw()
        #else:
            #print "%d -> %d: Not exists." % (val, new_val)
    s_timestep.on_changed(t_update)

    ### pwr_type selector
    pax = plt.axes([0.025, 0.40, 0.20, 0.35], axisbg=axcolor)
    radio1 = RadioButtons(pax, pwr_types, active=0)

    def type_update(val):
        t = get_t();
        g = get_group();

        pwr_type = str(val)
        set_pwr_type(val)

        print t, pwr_type, g

        y = result[t][pwr_type][g]
        l.set_ydata(y)
        plt.draw()

    radio1.on_clicked(type_update)

    ### group selector
    gax = plt.axes([0.025, 0.2, 0.15, 0.15], axisbg=axcolor)
    radio2 = RadioButtons(gax, groups, active=0)

    def g_update(val):
        t = get_t();
        pwr_type = get_pwr_type();

        g = int(val)
        set_group(val)

        print t, pwr_type, g

        y = result[t][pwr_type][g]
        l.set_ydata(y)
        plt.draw()
    radio2.on_clicked(g_update)

    #plt.show()
    '''

main()
