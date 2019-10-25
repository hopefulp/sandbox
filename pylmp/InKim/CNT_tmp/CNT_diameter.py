#!/home/noische/python

import sys
import math

usage = """
CNT_diameter.py "n m" prints the diameter.
CNT_diameter.py d prints the chiral index which satisfies the diameter.
"""

if len(sys.argv) < 2:
    print usage
    sys.exit(0)

parm = sys.argv[1:]
chiral_index = ""
for i in parm:
    chiral_index += str(i)
chiral_index = chiral_index.split()

if len(chiral_index) == 1:
    chiral_index = [float(s) for s in chiral_index]
    n = math.pi * chiral_index[0] / 0.246 / math.sqrt(3)
    print "- Chiral index of diameter %8.3f nm: %8.3f" % (chiral_index[0], n)
    sys.exit(0)

elif len(chiral_index) == 2:
    chiral_index = [int(s) for s in chiral_index]
    #chiral_index = [int(s) for s in chiral_index.split() if s.isdigit()]
    n, m = chiral_index
    d = 0.246 / math.pi * math.sqrt(n**2 + n*m + m**2)
    print "- Diameter of (%d, %d) CNT: %8.5f nm ( %8.5f in A )" % (chiral_index[0], chiral_index[1], d, d * 10)
    print "- Radius of (%d, %d) CNT: %8.5f nm ( %8.5f in A )" % (chiral_index[0], chiral_index[1], d/2, d * 5)
    sys.exit(0)

else:
    print "Wrong input for the chiral index (n, m): %s" % chiral_index
    sys.exit(0)
