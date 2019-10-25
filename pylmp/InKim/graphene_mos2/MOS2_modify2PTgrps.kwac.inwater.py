#!/home/noische/python

import sys
import os
import getopt
import time
import shutil
import math

from numpy import pi
import numpy as np
import bgf
import bgftools as bt
from grpfile import *
import nutils as nu

version = "160912"

def main(bgf_filename, grps_filename, margin = 0.0):

    # initialization
    grps_original_filename = grps_filename + ".original"

    # keep original grps file at the first run
    if not os.path.exists(grps_original_filename):
        shutil.copy(grps_filename, grps_original_filename)

    # read atoms from bgf file
    mybgf = bgf.BgfFile(bgf_filename)
    pbc = mybgf.CRYSTX[:3]
    swr = 3.270615945/2     # (wwr + swr) / 2 (defined by eq. from Pradeep Kumar's 2005 PRL)
    gwr = 3.057430885/2     # (wwr + cwr) / 2


    # wrap water molecules
    #mybgf = bt.periodicMoleculeSort(mybgf, pbc, selection="'WAT' in atom.rName", silent=False)
    #mybgf.saveBGF(bgf_filename + '.wrap')

    nu.warn("* Volume assignment from %s to %s" % (bgf_filename, grps_filename))

    '''
    Overall structure of MoS2
    -------------------------

    -----------Gra
    S_3b ----
    Mo      :rNo 2
    S_3a ----
            < layer_distance >
    S_3b ----
    Mo      :rNo 1
    S_3a ----
    -----------Gra
    '''

    # configure regions
    pbc_x = mybgf.CRYSTX[0]
    pbc_y = mybgf.CRYSTX[1]
    pbc_z = mybgf.CRYSTX[2]

    # group 1: MoS2
    avg_s3a_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 1")    # mos2 bottom
    avg_s3b_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 1")    # mos2 bottom
    avg_s3a_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3a' in atom.ffType and atom.rNo == 2")    # mos2 top
    avg_s3b_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'S_3b' in atom.ffType and atom.rNo == 2")    # mos2 top
    mos2_x = pbc_x
    mos2_y = pbc_y
    mos2_z = avg_s3b_z1 - avg_s3a_z1 + 2 * swr
    nu.warn("\t** MoS2 region **")
    nu.warn("\t\t- Effective MoS2 region dimensions: %8.3f %8.3f %8.3f" % (mos2_x, mos2_y, mos2_z))

    # group 2: graphene walls
    gwa_y = bt.atoms_average(mybgf, 'atom.y', selection="'GWA' in atom.rName")
    gwb_y = bt.atoms_average(mybgf, 'atom.y', selection="'GWB' in atom.rName")
    graphene_x = pbc_x
    graphene_y = gwb_y - gwa_y + 2 * gwr

    # group 3: inwater
    actual_distance = avg_s3a_z2 - avg_s3b_z1
    inwater_x = pbc_x
    inwater_y = graphene_y - 2 * margin
    inwater_z = actual_distance - 2 * swr    # effective interlayer distance
    graphene_z = pbc_z - mos2_z*2 - inwater_z   # graphene
    nu.warn("\t** In-water region **")
    nu.warn("\t\t- Actual interlayer distance: %8.3f" % actual_distance)
    nu.warn("\t\t- Effective water region dimensions: %8.3f %8.3f %8.3f" % (inwater_x, inwater_y, inwater_z))

    # group 4: reservoir
    outwater_x = pbc_x
    outwater_y = pbc_y - graphene_y
    outwater_z = pbc_z
    nu.warn("\t** Reservoir water region **")
    nu.warn("\t\t- dimensions: %8.3f %8.3f %8.3f" % (outwater_x, outwater_y, outwater_z))

    # find water locations
    all_molecules = bt.getMoleculeList(mybgf)
    all_ow = [atom.aNo for atom in mybgf.a if 'OW' in atom.ffType or 'HW' in atom.ffType]
    inwaters_ow = [atom.aNo for atom in mybgf.a if (atom.y > gwa_y + margin and atom.y < gwb_y - margin and 'OW' in atom.ffType)]
    inwaters_ow_nomargin = [atom.aNo for atom in mybgf.a if (atom.y > gwa_y and atom.y < gwb_y and 'OW' in atom.ffType)]

    if margin:
        nu.warn("\t- Margin specified: %8.3f" % margin)
        diff = len(inwaters_ow_nomargin) - len(inwaters_ow)
        if diff: nu.warn("\t- About %s water molecules are discarded out of %s atoms located near the border." % (diff, len(inwaters_ow_nomargin)))

    # record inwater aNo
    inwaters = []
    for ano in inwaters_ow:
        atom = mybgf.getAtom(ano)
        inwaters.append(atom.aNo)
        for i in atom.CONECT:
            inwaters.append(i)
        
    outwaters = [i for i in all_ow if not i in inwaters]
    nu.warn("** Distinguishing water positions **")
    nu.warn("\t- Found %d atoms in in-water region (%8.3f molecules)" % (len(inwaters), (len(inwaters)/3.0)))
    nu.warn("\t- Found %d atoms in reservoir region (%8.3f molecules)" % (len(outwaters), (len(outwaters)/3.0)))
    debug_selection = "same residue as (y > %8.3f and y < %8.3f and type OW)" % (gwa_y, gwb_y)
    nu.warn("VMD selection for inwaters: %s" % debug_selection)

    # calculate volume
    vol_mos2 = mos2_x * mos2_y * mos2_z
    vol_graphene = graphene_x * graphene_y * graphene_z
    vol_inwater = inwater_x * inwater_y * inwater_z
    vol_outwater = outwater_x * outwater_y * outwater_z

    # compute stats
    nu.warn("** Stats for confined water **")
    mass_inwater = 18.0154 * len(inwaters_ow)
    density_inwater = mass_inwater / vol_inwater / 6.022 * 10
    nu.warn("\t- Number: %d" % len(inwaters_ow))
    nu.warn("\t- Density: %8.5f" % density_inwater)

    # separate water group into inwater and outwater
    nu.warn("\tModifying grps file %s.." % sys.argv[0])
    g = grpfile(grps_original_filename)
    if len(g.grp) == 3:
        new_group_no = g.split_group(3, outwaters)
        g.grp[1]['volume'] = vol_mos2
        g.grp[2]['volume'] = vol_graphene
        g.grp[3]['volume'] = vol_inwater
        g.grp[4]['volume'] = vol_outwater
    elif len(g.grp) == 2:
        new_group_no = g.split_group(2, outwaters)
        g.grp[1]['volume'] = vol_mos2 + vol_graphene
        g.grp[2]['volume'] = vol_inwater
        g.grp[3]['volume'] = vol_outwater
    else:
        nu.die("Error on group numbers on %s" % grps_original_filename)
        
    #g.write(grps_filename, zip=False)
    g.write(grps_filename, zip=True)
        
    nu.warn("%s: Done." % sys.argv[0])
    ### end of code


if __name__ == "__main__":

    usage = """
MOS_modify2PTgrps_layer.py -b bgf_file -g grps_file

- This script modifies Group Atoms in *grps file
- This script should be run within Tod's LAMMPS-2PT pbs script.
"""

    if len(sys.argv) < 2:	
        print(usage);
        sys.exit(1)

    bgf_file = ""; grps_file = ""; margin = 0.0;

    options, args = getopt.getopt(sys.argv[1:], 'hb:g:m:', ['help', 'bgffile=', 'grpsfile=', 'margin='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            sys.exit(0);
        elif option in ('-b', '--bgffile'):
            bgf_file = value
        elif option in ('-g', '--grpsfile'):
            grps_file = value
        elif option in ('-m', '--margin'):
            margin = float(value)

    main(bgf_file, grps_file, margin=margin)
