###### keywords to deal with file

### Molden format from Q-Chem to draw using molden
qcout_molden2="MOLDEN-FORMATTED"        # 2nd (last) block of MOLDEN "MOLDEN-FORMATTED INPUT FILE
qcout_molden1="MOLDEN FORMAT"           # 1st molden block; start by    "MOLECULAR ORBITALS IN MOLDEN FORMAT"
qcout_molden1b="MOLDEN-FORMAT"          # 1st molden bolck: end by      "MOLDEN-FORMAT MOLECULAR ORBITALS"

### extract certain bands from BAND.dat from vaspkit

def get_kw(jobtype):
    if jobtype == 'molden':
        kw1 = qcout_molden1
        kw2 = qcout_molden1b
    elif jobtype == 'band':
        kw1 = 'Band-Index'
        kw2 = None
    return kw1, kw2


