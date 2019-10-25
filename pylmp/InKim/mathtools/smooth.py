import numpy as np

def smooth(x,window_len=11,window='hanning'):
    """from http://scipy.github.io/old-wiki/pages/Cookbook/SignalSmooth
    """ 
     
    if x.ndim != 1:
        print("smooth only accepts 1 dimension arrays.")
        raise ValueError

    if x.size < window_len:
        print("Input vector needs to be bigger than window size.")
        raise ValueError
        

    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
        raise ValueError
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')

    return y # original: WARNING: len(input) != len(output)
    #print(np.ceil(window_len/2-1))
    #print(-np.ceil(window_len/2))
    #print(len(y[np.ceil(window_len/2-1):-np.ceil(window_len/2)]))
    #return y[int(np.ceil(window_len/2-1)):int(-np.ceil(window_len/2)+1)]

