SYSTEM = zzrpx2-afm-rlx-re

MAGMOM = 1.5 1.5  999*0
# continuation
ISTART = 0
ICHARG = 0
ISPIN = 2

# precision-1
ENCUT = 500
PREC = high
ISMEAR = 0 ; SIGMA = 0.05
NELMIN = 4 ; NELM = 500
EDIFF = 1E-5 ; EDIFFG = -0.025

# precision-2
GGA_COMPAT=.FALSE.
#VOSKOWN = 1
ADDGRID = .TRUE.

# mixing
LMAXMIX = 4
#IMIX = 4
#AMIX = 0.2
#BMIX = 0.0001
#AMIX_MAG = 0.8
#BMIX_MAG = 0.0001

# performance
ALGO = Fast
#IALGO=48
NSIM = 4
NPAR = 4
LREAL = Auto
LPLANE = .TRUE.
LSCALAPACK = .FALSE.

# geo-opt
#NSW = 999 ; ISIF = 3
#IBRION = 2
#POTIM = 0.3

# aimd - nvt
#NSW = 5000 ; ISIF = 0
#IBRION = 0
#POTIM = 2.0
#TEIN = 300
#TEBEG = 300
#TEEND = 300
#SMASS = 0.05
#ISYM = 0

# log
LORBIT = 11
LAECHG=.TRUE.
LWAVE = .TRUE.
LCHARG = .TRUE.

# functional
GGA = RP

# U-correction
LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = 2 -1 -1 -1
LDAUU = 6.5 0 0 0
LDAUJ = 1 0 0 0

# D2 correction
#LVDW = .TRUE.

# hse06
#GGA = PE
#LHFCALC = .TRUE.
#ALGO = D; TIME = 0.4
#HFSCREEN = 0.207
#PRECFOCK = F
#NKRED = 2
#LTHOMAS = .TRUE.
#Block NPAR tag

# pbe0
#GGA = PE
#LHFCALC = .TRUE.
#ALGO = D; TIME = 0.4
#ENCUTFOCK = 0
#NKRED = 2
#Block NPAR tag

# b3lyp
#GGA = B3
#LHFCALC = .TRUE.
#ALGO = D; TIME = 0.5
#AEXX = 0.2
#AGGAX = 0.72
#AGGAC = 0.81
#ALDAC = 0.19
#Block NPAR tag

# dos
EMIN = -60
EMAX = 20
NEDOS = 8000

# orbital
#NBANDS = 112
#LPARD = .TRUE.
#IBAND = 91 92 93 94 95 96 97 98
#KPUSE = 1 2 3 4 5 6
#LSEPB = .TRUE.
#LSEPK = .TRUE.
