import numpy as np

'''
VASP analysis tool
    calculate band charge amount from EIGENVAL
'''

def defect_band_electron_count(eigenval_file, target_band):
    """
    Compute total electron occupation of a specific band
    from EIGENVAL (ISPIN = 1).
    
    Parameters
    ----------
    eigenval_file : str
        Path to EIGENVAL
    target_band : int
        Band index (1-based, as printed in EIGENVAL)

    Returns
    -------
    float
        Total electrons occupying that band
    """

    with open(eigenval_file, 'r') as f:
        lines = f.readlines()

    # Line 6 contains: nelect, nkpoints, nbands
    header = lines[5].split()
    nkpoints = int(header[1])
    nbands = int(header[2])

    idx = 7  # start of first k-point block (standard format)
    total_electrons = 0.0

    for k in range(nkpoints):
        # Read k-point line
        k_line = lines[idx].split()
        weight = float(k_line[3])
        idx += 1

        # Loop over bands
        for b in range(1, nbands + 1):
            band_line = lines[idx].split()
            band_index = int(band_line[0])
            energy = float(band_line[1])
            occ = float(band_line[2])

            if band_index == target_band:
                total_electrons += 2.0 * occ * weight
                print(target_band, occ, weight)

            idx += 1

        idx += 1  # skip blank line between k-points


    return total_electrons
