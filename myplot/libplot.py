'''
lib for plot
    extract x, y
'''

def extract_xy_nfile(files, ix, iy):

    xnc = []
    ync = []
    for f1 in files:
        x_vals = []
        y_vals = []
        with open(f1, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    x, y = map(float, line.split())
                    x_vals.append(x)
                    y_vals.append(y)

        xnc.append(x_vals)
        ync.append(y_vals)
    return xnc, ync