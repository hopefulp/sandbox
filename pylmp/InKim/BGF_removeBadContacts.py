#!/home/noische/Enthought/Canopy_64bit/User/bin/python

# Python Modules
import sys
import re
import getopt
import optparse

# BGF modules
import bgftools

__version__ = "20160405"

def remove_bad_contacts(bgf_file, out_file='', thresh=2.0, silent=False):

    bgftools.remove_bad_contacts(bgf_file, out_file=out_file, thresh=thresh, silent=False)

if __name__ == "__main__":

    # Globals
    option = ""; args = ""; bgf_file = ""; out_file = ""; thresh = "";
    usage = """
    removeBadContacts.py: read coordinate data from the LAMMPS trajectory file
               write the data to the original BGF file
    Usage: removeBadContacts.py -b bgf_file -t distance_threshold -o out_file
    """

    print("Requested options: %s" % sys.argv)

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:t:o:', ['help','bgf=','thresh=','out='])
    print("Requested options: " + str(options))
    for option, value in options:
        if option in ('-h', '--help'):
            print usage
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = value
        elif option in ('-t', '--thresh'):
            thresh = float(value)
        elif option in ('-o', '--out'):
            out_file = value

    if not thresh:
        thresh = 2.4

    # main call
    remove_bad_contacts(bgf_file, out_file, thresh)

