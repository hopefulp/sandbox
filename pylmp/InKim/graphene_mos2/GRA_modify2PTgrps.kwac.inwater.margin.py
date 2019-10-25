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

def main(bgf_filename, grps_filename, margin):

    # initialization
    grps_original_filename = grps_filename + ".original"
    warn = ""

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

    warn += "\n* Volume assignment from %s to %s\n" % (bgf_filename, grps_filename)

    '''
    Overall structure of MoS2
    -------------------------

    -----------GRB
            < layer_distance >
    -----------GRA
    '''

    # configure regions
    pbc_x = mybgf.CRYSTX[0]
    pbc_y = mybgf.CRYSTX[1]
    pbc_z = mybgf.CRYSTX[2]

    # group 1: Graphene interlayer
    avg_z1 = bt.atoms_average(mybgf, 'atom.z', selection="'C_R' in atom.ffType and 'GRA' in atom.rName")    # gra bottom
    avg_z2 = bt.atoms_average(mybgf, 'atom.z', selection="'C_R' in atom.ffType and 'GRB' in atom.rName")    # gra top
    gra_x = pbc_x
    gra_y = pbc_y
    gra_z = 2 * gwr
    warn += "\t** Graphene interlayer region **\n"
    warn += "\t\t- Effective graphene region dimensions: %8.3f %8.3f %8.3f\n" % (gra_x, gra_y, gra_z)

    # group 2: graphene walls
    gwa_y = bt.atoms_average(mybgf, 'atom.y', selection="'GWA' in atom.rName")
    gwb_y = bt.atoms_average(mybgf, 'atom.y', selection="'GWB' in atom.rName")
    wall_x = pbc_x
    wall_y = gwb_y - gwa_y + 2 * gwr

    # group 3: "real" inwater (center inwater)
    actual_distance = avg_z2 - avg_z1
    inwater_x = pbc_x
    inwater_y = wall_y - 2 * margin
    inwater_z = actual_distance - 2 * gwr    # effective interlayer distance
    wall_z = pbc_z - gra_z*2 - inwater_z   # graphene
    warn += "\t** In-water region **\n"
    warn += "\t\t- Actual interlayer distance: %8.3f\n" % actual_distance
    warn += "\t\t- Effective water region dimensions: %8.3f %8.3f %8.3f\n" % (inwater_x, inwater_y, inwater_z)

    # group 4: marginal area of inwater (entrance)
    marginal_inwater_x = pbc_x
    marginal_inwater_y = 2 * margin
    marginal_inwater_z = inwater_z

    # group 5: marginal area of reservoir water (entrance)
    marginal_outwater_x = pbc_x
    marginal_outwater_y = 2 * margin
    marginal_outwater_z = pbc_z

    # group 6: "real" reservoir (waterly region)
    outwater_x = pbc_x
    outwater_y = pbc_y - wall_y - 2 * margin
    outwater_z = pbc_z
    warn += "\t** Reservoir water region **\n"
    warn += "\t\t- dimensions: %8.3f %8.3f %8.3f\n" % (outwater_x, outwater_y, outwater_z)

    # find water locations
    #all_molecules = bt.getMoleculeList(mybgf)
    all_ow = [atom.aNo for atom in mybgf.a if 'OW' in atom.ffType]
    inwaters_ow = [atom.aNo for atom in mybgf.a if (atom.y > gwa_y + margin and atom.y < gwb_y - margin and 'OW' in atom.ffType)]   # real inwater
    marginal_inwaters_ow = [atom.aNo for atom in mybgf.a if ((atom.y > gwa_y and atom.y < gwa_y + margin) or (atom.y > gwb_y - margin and atom.y < gwb_y)) and 'OW' in atom.ffType] # marginal inwater
    inwaters_ow_nomargin = [atom.aNo for atom in mybgf.a if (atom.y > gwa_y and atom.y < gwb_y and 'OW' in atom.ffType)]    # inwater = real + marginal
    marginal_outwaters_ow = [atom.aNo for atom in mybgf.a if ((atom.y > gwa_y - margin and atom.y < gwa_y) or (atom.y > gwb_y and atom.y < gwb_y + margin)) and 'OW' in atom.ffType]    # marginal outwater
    outwaters_ow = [atom.aNo for atom in mybgf.a if not atom.aNo in inwaters_ow and not atom.aNo in marginal_outwaters_ow and not atom.aNo in marginal_inwaters_ow and 'OW' in atom.ffType]    # real outwater

    print("inwaters_vmd:  same residue as (y > {gwa_y} + {margin} and y < {gwb_y} - {margin} and type OW)".format(**vars()))
    print("marginal_inwaters_vmd:  same residue as (((y > {gwa_y} and y < {gwa_y} + {margin}) or (y > {gwb_y} - {margin} and y < {gwb_y})) and type OW)".format(**vars()))
    print("inwaters_vmd_nomargin:  same residue as (y > {gwa_y} and y < {gwb_y} and type OW)".format(**vars()))
    print("marginal_outwaters_vmd:  same residue as (((y > {gwa_y} - {margin} and y < {gwa_y}) or (y > {gwb_y} and y < {gwb_y} + {margin})) and type OW)".format(**vars()))

    if len(inwaters_ow) + len(marginal_inwaters_ow) != len(inwaters_ow_nomargin):
        warn += "Number of water molecules mismatch in inwater!\n"

    warn += "\t- Margin specified: %8.3f\n" % margin

    # record aNo
    def get_atoms(lst):
        result = []
        for ano in lst:
            atom = mybgf.getAtom(ano)
            result.append(atom.aNo)
            for i in atom.CONECT:
                result.append(i)
        return result

    inwaters = get_atoms(inwaters_ow)
    marginal_inwaters = get_atoms(marginal_inwaters_ow)
    marginal_outwaters = get_atoms(marginal_outwaters_ow)
    outwaters = get_atoms(outwaters_ow)
        
    if not (len(inwaters) % 3 == 0 and len(marginal_inwaters) % 3 == 0 and len(marginal_outwaters) % 3 == 0 and len(outwaters) % 3 == 0):
        warn += "----- Suspicious water molecules division found!!!! -----\n"

    # calculate volume
    vol_gra = gra_x * gra_y * gra_z
    vol_graphene = wall_x * wall_y * wall_z
    vol_inwater = inwater_x * inwater_y * inwater_z
    vol_marginal_inwater = marginal_inwater_x * marginal_inwater_y * marginal_inwater_z
    vol_outwater = outwater_x * outwater_y * outwater_z
    vol_marginal_outwater = marginal_outwater_x * marginal_outwater_y * marginal_outwater_z

    vol_total = pbc_x * pbc_y * pbc_z
    vol_sum = vol_gra + vol_graphene + vol_inwater + vol_marginal_inwater + vol_outwater + vol_marginal_outwater

    # compute stats
    warn += "\n\n\t**** Stats for confined water ****\n"
    warn += "\t- Number:\n"
    warn += "\t\t real inwater / marginal inwater / marginal outwater / real outwater: %d %d %d %d\n" % (len(inwaters_ow), len(marginal_inwaters_ow), len(marginal_outwaters_ow), len(outwaters_ow))
    warn += "\t\t (%5.1f %%  %5.1f %%  %5.1f %%  %5.1f %%)\n" % (float(len(inwaters_ow))/len(all_ow)*100, float(len(marginal_inwaters_ow))/len(all_ow)*100, float(len(marginal_outwaters_ow))/len(all_ow)*100, float(len(outwaters_ow))/len(all_ow)*100)

    warn += "\t- Volumes: \n"
    warn += "\t\tvol_gra: %.3f = %.3f * %.3f * %.3f\n" % (vol_gra, gra_x, gra_y, gra_z)
    warn += "\t\tvol_graphene: %.3f = %.3f * %.3f * %.3f\n" % (vol_graphene, wall_x, wall_y, wall_z)
    warn += "\t\tvol_inwater: %.3f = %.3f * %.3f * %.3f\n" % (vol_inwater, inwater_x, inwater_y, inwater_z)
    warn += "\t\tvol_marginal_inwater: %.3f = %.3f * %.3f * %.3f\n" % (vol_marginal_inwater, marginal_inwater_x, marginal_inwater_y, marginal_inwater_z)
    warn += "\t\tvol_marginal_outwater: %.3f = %.3f * %.3f * %.3f\n" % (vol_marginal_outwater, marginal_outwater_x, marginal_outwater_y, marginal_outwater_z)
    warn += "\t\tvol_outwater: %.3f = %.3f * %.3f * %.3f\n" % (vol_outwater, outwater_x, outwater_y, outwater_z)

    if vol_sum != vol_total:
        warn += "\n\t\t----- Suspicious volume division found!!!! Density might be different -----\n"
        warn += "\t\tTotal volume: %8.3f, Sum of divisions: %8.3f\n" % (vol_total, vol_sum)

    warn += "\t- Density:\n"
    density_inwater = 18.0154 * len(inwaters_ow) / vol_inwater / 6.022 * 10
    density_marginal_inwater = 18.0154 * len(marginal_inwaters_ow) / vol_marginal_inwater / 6.022 * 10
    density_marginal_outwater = 18.0154 * len(marginal_outwaters_ow) / vol_marginal_outwater / 6.022 * 10
    density_outwater = 18.0154 * len(outwaters_ow) / vol_outwater / 6.022 * 10
    warn += "\t\t real inwater / marginal inwater / marginal outwater / real outwater: %.3f %.3f %.3f %.3f\n" % (density_inwater, density_marginal_inwater, density_marginal_outwater, density_outwater)

    # separate water group into inwater and outwater
    warn += "Modifying grps file %s..\n" % sys.argv[0]
    g = grpfile(grps_original_filename) # suppose there are two groups in the grps file
    marginal_inwater_group_no = g.split_group(2, marginal_inwaters)
    marginal_outwater_group_no = g.split_group(2, marginal_outwaters)
    outwater_group_no = g.split_group(2, outwaters)
    g.grp[1]['volume'] = vol_gra + vol_graphene
    g.grp[2]['volume'] = vol_inwater
    g.grp[marginal_inwater_group_no]['volume'] = vol_marginal_inwater
    g.grp[marginal_outwater_group_no]['volume'] = vol_marginal_outwater
    g.grp[outwater_group_no]['volume'] = vol_outwater
        
    #g.write(grps_filename, zip=False)
    g.write(grps_filename, zip=True)
        
    warn += "%s: Done." % sys.argv[0]

    #nu.warn(warn)
    print("Numbers: real inwater / marginal inwater / marginal outwater / real outwater: %d %d %d %d" % (len(inwaters_ow), len(marginal_inwaters_ow), len(marginal_outwaters_ow), len(outwaters_ow)))
    print("Density: real inwater / marginal inwater / marginal outwater / real outwater: %.3f %.3f %.3f %.3f" % (density_inwater, density_marginal_inwater, density_marginal_outwater, density_outwater))
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

    main(bgf_file, grps_file, margin)
