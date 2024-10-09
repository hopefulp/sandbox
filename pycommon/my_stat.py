#### statistics 

import numpy as np
from ase.calculators.calculator import Parameters

def stat_2col(x, y, scale=1):
    p = Parameters()
    x = np.array(x)
    y = np.array(y)

    diff = x - y
    #for i, (f_true, f_hypo) in enumerate(zip(x, y)):
    #    print(f" {i}: true: {f_true}, hypo: {f_hypo}")
    maxres  = np.max(np.absolute(diff))
    imaxres = np.argmax(diff)
    meanx   = np.mean(x)
    meany   = np.mean(y)
    bias_sq = (meanx - meany)**2                # Bias(X-Y)**2

    varx    = np.var(x)                         # Variance(X)
    vary    = np.var(y)                         # Variance(Y)
    cov     = np.mean((x-meanx) * (y-meany))    # Covariance
    varx_y  = varx + vary - 2 * cov             # variance(X-Y)
    mse     = bias_sq + varx_y                  # MSE
    rmse    = np.sqrt(mse)
    r_pearson = cov/np.sqrt(varx*vary)

    if scale != 1:
        mse     *= (scale**2)
        bias_sq *= (scale**2)
        varx_y  *= (scale**2)
        varx    *= (scale**2)
        vary    *= (scale**2)
        cov     *= (scale**2)
        rmse    *= scale
        maxres  *= scale
    p.mse       = mse
    p.bias_sq   = bias_sq
    p.varx_y    = varx_y
    p.varx      = varx
    p.vary      = vary
    p.cov       = cov
    p.r_corr    = r_pearson
    p.rmse      = rmse
    p.maxres    = maxres
    p.imaxres   = imaxres

    return p
    #return mse, bias_sq, varx_y, r_pearson, varx, vary, cov

def get_MBD_1D(loc=0, scale=1, size = None):
    """Generates size samples of maxwell
    loc = mean
    scale = sigma
    normal distribution for v_x, v_y, v_z makes MB distribition of v
    """
    vx = np.random.normal(loc=loc, scale=scale, size=size)
    vy = np.random.normal(loc=loc, scale=scale, size=size)
    vz = np.random.normal(loc=loc, scale=scale, size=size)
    return vx, vy, vz

def get_MBD_v(loc=0, scale=1, size = None):
    vx, vy, vz = get_MBD_1D(loc=0, scale=1, size = None)
    return np.sqrt(vx*vx + vy*vy + vz*vz)
 
