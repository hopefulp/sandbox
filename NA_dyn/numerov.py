import numpy as np

def numerov(xgrid, y_0, y_1, k, S):
    """
    Inputs
        xgrid : x (1-D array) in real scale with reduced dimension
        y_0   : y, initial condition 1
    Output
        y : solution of the differential equation (1-D array)
    """
    # initialize y
    ngrid = len(xgrid)
    h = xgrid[1] - xgrid[0]
    y = np.zeros(ngrid)
    y[0] = y_0
    y[1] = y_1

    # main loop: evaluate y[j]
    for j in np.arange(2, ngrid):

        y1 = y[j-2]; y2 = y[j-1]
        k1 = k[j-2]; k2 = k[j-1]; k3 = k[j]
        s1 = S[j-2]; s2 = S[j-1]; s3 = S[j]

        term_S = 1/12. * h**2 * (s3 + 10*s2 + s1)
        term_3 =      (1 + 1/12. * h**2 * k3**2)
        term_2 = -2 * (1 - 5/12. * h**2 * k2**2) * y2
        term_1 =      (1 + 1/12. * h**2 * k1**2) * y1

        y3 = (term_S - term_2 - term_1) / term_3
        y[j] = y3

    return y

def numerov_k2(xgrid, y_0, y_1, k, S):
    """
    To deal with the sign of (+-)k^2, k should be k(x)^2 in the input
    Inputs
        xgrid : x (1-D array) in real scale with reduced dimension
        y_0   : y, initial condition 1
    Output
        y : solution of the differential equation (1-D array)
    """
    # initialize y
    ngrid = len(xgrid)
    h = xgrid[1] - xgrid[0]
    y = np.zeros(ngrid)
    y[0] = y_0
    y[1] = y_1

    # main loop: evaluate y[j]
    for j in np.arange(2, ngrid):

        y1 = y[j-2]; y2 = y[j-1]
        k1 = k[j-2]; k2 = k[j-1]; k3 = k[j]
        s1 = S[j-2]; s2 = S[j-1]; s3 = S[j]

        term_S = 1/12. * h**2 * (s3 + 10*s2 + s1)
        term_3 =      (1 + 1/12. * h**2 * k3)
        term_2 = -2 * (1 - 5/12. * h**2 * k2) * y2
        term_1 =      (1 + 1/12. * h**2 * k1) * y1

        y3 = (term_S - term_2 - term_1) / term_3
        y[j] = y3

    return y

def numerov_inward(xgrid, y_0, y_1, k, S):

    # initialize y
    ngrid = len(xgrid)
    h = np.abs(xgrid[1] - xgrid[0])
    y = np.zeros(ngrid)
    y[-1] = y_0
    y[-2] = y_1

    # main loop: evaluate y[j]
    for j in np.arange(2, ngrid):
        #print (j, 1-j, 0-j, -1-j)

        y2 = y[-j];   y3 = y[-j+1]
        k1 = k[-j-1]; k2 = k[-j]; k3 = k[-j+1]
        s1 = S[-j-1]; s2 = S[-j]; s3 = S[-j+1]
        term_S = 1/12. * h**2 * (s3 + 10*s2 + s1)
        term_3 =      (1 + 1/12. *   h**2 * k3**2) * y3
        term_2 = -2 * (1 - 5/12. *   h**2 * k2**2) * y2
        term_1 =      (1 + 1/12. *   h**2 * k1**2)

        y1 = (term_S - term_2 - term_3) / term_1
        y[-j-1] = y1

    return y


