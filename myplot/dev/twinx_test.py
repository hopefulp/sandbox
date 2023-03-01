"""
    2d plot library
    draw_2d
    fdraw   where f == "file"
    barplot
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import sys
import my_chem 
import numpy as np
from common import *

#plt.switch_backend('agg')
size_title = 20
size_label = 18
#size_tick = 15
text_x = 0.75
text_y = 0.8
text_twinx_x = 0.8
text_twinx_y = 0.9

### functions
#   my_font
#   common_figure: will be moved to myplot_default and deleted
#   common_figure_after: needed after plt.plot()
#   draw_dots_two: option: twinx
#   xtitle_font
#   mplot_twinx
#   mplot_nvector
#   draw_histogram
#

#from myplot_default import *

def mplot_twinx(x, y, iy_right, title=None, xlabel=None, ylabel='E [eV]', legend=None, Lsave=False, Colors=None):
    '''
    called from "amp_plot_stat.py"
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    len(ylabel) == 2
    '''
    #fig, ax = common_figure(ncolor = len(y))
    fig, ax = plt.subplots(figsize=(12,8))
    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    ### function: mplot_twinx
    plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=25)
    if isinstance(ylabel, str): ylabel1 = ylabel2 = ylabel
    elif isinstance(ylabel, list):
        ylabel1 = ylabel[0]
        ylabel2 = ylabel[1]
    plt.ylabel(ylabel1, fontsize=25, color='r')
    ax.tick_params(axis='y', colors='r')
    #ax.xaxis.set_major_locator(plt.NullLocator())
    print(f"x, y shape:: {np.array(x).shape} {np.array(y).shape} and ylabel {ylabel} in {whereami()}")
    ax2=ax.twinx()
    ax2.set_ylabel(ylabel2, fontsize=25, color='g')
    #ax2.set_ylim(0,0.2)
    pls=[]
    print(f"{iy_right} {len(ys)} {whereami()}")
    for i in range(len(ys)):
        if i not in iy_right:
            print(f"{i} draw on the left")
            plt.yticks(color='r')
            p1, = ax.plot(x, ys[i,:], '-', color='r',  label=legend[i])
            pls.append(p1)
        else:
            print(f"{i} draw on the right")
            plt.yticks(color='g')
            p2, = ax2.plot(x, ys[i,:], '-', color='g', label=legend[i])
            pls.append(p2)
    ax.legend(loc=3)            # 2
    ax2.legend(loc=4)           # 1
    plt.legend()
    plt.show()
    return 0

