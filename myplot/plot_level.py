"""
    2d plot library
    draw_2d
    fdraw   where f == "file"
    barplot
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from cycler import cycler
import sys
import my_chem 
import numpy as np
from common import *
#plt.switch_backend('agg')
size_title = 20
size_label = 18
size_tick = 15
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

def my_font(pack='amp'):
    if pack == 'amp':
        font={  'family': 'normal',
                'weight': 'bold',
                'size'  : 22,
                }
    else:
        print("package: %s is not included" % pack)
        sys.exit(0)
    mpl.rc('font', **font)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    return
### choice between "import myplot_default"|call common_figure
def common_figure():
    fig = plt.figure(figsize=(15,10))
    ax = plt.axes()
    mpl.rcParams.update({'font.size':25})
    ax.tick_params(axis='both', which='major', labelsize=25)
    #ax.tick_params(axis='x', labelsize=30)

    #custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))       # Figure 8(a)
    custom_cycler = (cycler(color=['r','g']))                                          # Figure 8(b)
    #custom_cycler = (cycler(color=['darkcyan','b']))                                   # Figure S16(a)
    #custom_cycler = (cycler(color=['r','darkcyan']))                                    # Figure S16(b)
    #custom_cycler = (cycler(color=['orange','m','g','b']))
    ax.set_prop_cycle(custom_cycler)
    return fig, ax

### draw_dots_two was upgraded to twinx
def draw_2subdots(y, h, title, suptitle, Ltwinx=None, escale=1.0,Colors=['r','b','o'], Ldiff=True):
    '''
    this makes error in serial plotting
    '''
    fig, ax = common_figure()
    escale = 1.0
    nlen = len(y)
    h_conv = np.array(h) * escale        # escale = my_chem.ev2kj
    y_conv = np.array(y) * escale
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())
    max_res = abs(max(diff, key=abs))
    #max_res = max(diff, key=abs)
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text_pos_x = nlen*0.85                  # 0.85, 0.2
    text_pos_y = max(y_conv)*0.2
    text="E_rms(test) = {:7.3f}\nE_maxres = {:7.3f}".format(rmse, max_res)

    if Colors:  color = Colors.pop(0)       #'tab:' + Colors.pop(0)
    else:       color='tab:green'
    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    #mpl.rcParams.update({'font.size':22})
    plt.title(title, fontsize=20)
#    plt.suptitle(suptitle, x=0.5, y=0.92, va='top', fontsize=18)
    plt.suptitle(suptitle, fontsize=10)
    if escale == 1.0:
        ax.set_ylabel('PE(eV)', color='b', fontsize=15)
    elif escale == my_chem.ev2kj:
        ax.set_ylabel('PE(kJ/mol)', fontsize=15)
    plt.xlabel('data', fontsize=15)
    #ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='y', labelcolor='b', labelsize=10)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    if Ltwinx:
        ax2=ax.twinx()
        ax2.set_ylabel("Difference(eV)", color='g')
        #ax2.plot(x, ys[1,:], 'o-', color=color)
        ax2.tick_params(axis='y', labelcolor='g', labelsize=10)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    p1  = ax.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = ax.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    print(y_conv)
    
    if Ltwinx:
        if Ldiff:
            p3, = ax2.plot(range(nlen), diff, c='g', label='difference')
            plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
        else:
            plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_twinx_x, text_twinx_y, text, fontsize=10, transform=ax.transAxes)
    else:
        #p3, = plt.plot(range(nlen), diff, c='g', label='difference')
        #plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
        plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_x, text_y, text, fontsize=10, transform=ax.transAxes)
    plt.plot(range(nlen), ones)

    plt.show()
    return rmse, max_res
### draw_dots_two was upgraded to twinx
### still used by amp_run.py
def draw_amp_twinx(y, h, title, suptitle, natom=1, Ltwinx=None, escale=1.0,Colors=['r','b','o'], Ldiff=True):
    '''
    this makes error in serial plotting
    '''
    fig, ax = common_figure()
    escale = 1.0
    nlen = len(y)
    h_conv = np.array(h) * escale        # escale = my_chem.ev2kj
    y_conv = np.array(y) * escale
    ymin = min(y_conv)
    ymax = max(y_conv)
    y_width = ymax - ymin
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())/natom                      # eV/atom
    max_res = abs(max(diff, key=abs))/natom                     # eV/atom
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text_pos_x = nlen*0.85                  # 0.85, 0.2
    text_pos_y = max(y_conv)*0.2
    text="E_rms(test) = {:7.3f} eV/atom\nE_maxres   = {:7.3f} eV/atom".format(rmse, max_res)

    if Colors:  color = Colors.pop(0)       #'tab:' + Colors.pop(0)
    else:       color='tab:green'
    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    #mpl.rcParams.update({'font.size':22})
    plt.title(title+'\n', fontsize=size_title)
    #suptitle=suptitle+'\n'
    plt.suptitle(suptitle, x=0.5, y=0.96, va='top', fontsize=18)
    #plt.suptitle(suptitle, fontsize=size_tick)
    if escale == 1.0:
        ax.set_ylabel('PE(eV)', color='b', fontsize=size_label)
        ax.set_ylim(ymin-1, ymax+0.25)
    elif escale == my_chem.ev2kj:
        ax.set_ylabel('PE(kJ/mol)', fontsize=size_tick)
    plt.xlabel('data', fontsize=size_label)
    #ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='y', labelcolor='b', labelsize=size_tick)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    if Ltwinx:
        ax2=ax.twinx()
        ax2.set_ylabel("Difference(eV)", color='g')
        ax2.set_ylim(-0.001, 0.01)
        #ax2.plot(x, ys[1,:], 'o-', color=color)
        ax2.tick_params(axis='y', labelcolor='g', labelsize=size_tick)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    #print(y_conv, h_conv)
    p1  = ax.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = ax.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    if Ltwinx:
        if Ldiff:
            p3, = ax2.plot(range(nlen), diff, c='g', label='difference')
            plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.2))
        else:
            plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_twinx_x, text_twinx_y, text, fontsize=size_tick, transform=ax.transAxes)
    else:
        #p3, = plt.plot(range(nlen), diff, c='g', label='difference')
        #plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
        plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_x, text_y, text, fontsize=10, transform=ax.transAxes)
    plt.plot(range(nlen), ones)

    plt.show()
    return rmse, max_res





def xtitle_font(tit):
    st = "\'{}\', fontsize=20".format(tit)
    print(st)
    return st


### used for md.ene, normal data file
#def mplot_twinx(x, y, dx=1.0, Title=None, Xtitle=None, Ytitle=None, Ylabels=None, Lsave=False, Colors=None):
def mplot_twinx(x, y, iy_right, title=None, xlabel=None, ylabel=None, legend=None, Lsave=False, Colors=None):
    '''
    called from "amp_plot_stat.py"
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    len(ylabel) == 2
    '''
    fig, ax = common_figure()
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
    #print(f"x, y shape:: {np.array(x).shape} {np.array(y).shape} and ylabel {ylabel} in {whereami()}")
    ax2=ax.twinx()
    ax2.set_ylabel(ylabel2, color='g')
    #ax2.set_ylim(0,0.2)
    pls=[]
    print(iy_right, len(ys))
    for i in range(len(ys)):
        if i in iy_right: 
            #if Colors:  color = Colors.pop(i)       #'tab:' + Colors.pop(0)
            #else:       color='tab:green'
            plt.yticks(color='g')
            p2, = ax2.plot(x, ys[i,:], 'x-', color='g', label=legend[i])
            pls.append(p2)
        else:
            #ax2.tick_params(axis='y')
            #if Colors:  color = Colors.pop(i)       #color = 'tab:' + Colors.pop(0)
            #else:       color = 'tab:red'
            #print(f"shape of x, ys[i] = {np.array(x).shape} {ys[i,:].shape}")
            plt.yticks(color='r')
            p1, = ax.plot(x, ys[i,:], 'o-', color='r',  label=legend[i])
            pls.append(p1)
    #plt.legend(pls, Ylabels, loc=2)
    #ax.legend(loc=3)            # 2
    #ax2.legend(loc=4)           # 1
    #plt.legend()
    #common_figure_after()
    #x_ticks = ['PP', 'PPP', 'PNP', 'PNP-bridged']
    #plt.xticks(x_ticks)
    #ax.set_xticklabels(x_ticks)
    plt.locator_params(axis='x', nbins=10)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0


def mplot_levels(x, y, dx=1.0, title=None, xlabel=None, ylabel=None, legend=None,Lsave=False, Colors=None):
    '''
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    '''
    fig, ax = common_figure()
    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    print(f"x={len(x)} y={ys.shape}")
    plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    if ys.ndim == 1:
        plt.plot(x, y, 'bo-')
        #plt.scatter(x, y)
    elif ys.ndim == 2:
        for i in range(len(ys)):
            plt.plot(x,ys[i,:], 'o-', label=legend[i] )
    else:
        print(f"Error:: obscure in y-dim {ys.ndim}")
    ### ADD LEGEND
    #plt.legend(loc=2)                # locate after plot
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    #ax.yaxis.set_major_locator(ticker.MultipleLocator())
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0



def mplot_nvector(x, y, dx=1.0, title=None, xlabel=None, ylabel=None, legend=None,Lsave=False, Colors=None):
    '''
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    '''
    fig, ax = common_figure()
    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    print(f"x={len(x)} y={ys.shape}")
    plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    if ys.ndim == 1:
        plt.plot(x, y, 'bo-')
        #plt.scatter(x, y)
    elif ys.ndim == 2:
        for i in range(len(ys)):
            plt.plot(x,ys[i,:], 'o-', label=legend[i] )
    else:
        print(f"Error:: obscure in y-dim {ys.ndim}")
    ### ADD LEGEND
    #plt.legend(loc=2)                # locate after plot
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    #ax.yaxis.set_major_locator(ticker.MultipleLocator())
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0






def mplot_nvector_v1(x, y, dx=1.0, Title=None, Xtitle=None, Ytitle=None, Ylabels=None, Lsave=False, Colors=None, Ltwinx=None):
    '''
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[multi],[multi],...size]
    '''
    if Ltwinx:
        ax2 = ax.twinx()
        ax2.set_ylabel("Kinetic energy(kJ/mol)")
        ax2.tick_params(axis='y', labelcolor='g', labelsize=10)

    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    #print("hmm: {}".format(xsize))
    print(f"{x} :: {ys}")
    #plt.xticks(np.arange(min(x), max(x)+1, int(max(x)/dx)))
    #plt.xticks(np.arange(min(x), max(x)+1))
    #if tag=='x-sub':
    #    #xlabels = [item.get_text() for item in ax.get_xticklabels()]
    #    xlabels = tag_value
    #    ax.set_xticklabels(xlabels)

    plt.title(Title)
    if Xtitle:
        plt.xlabel(Xtitle, fontsize=15)
    plt.ylabel(Ytitle, fontsize=15)
    #ax.xaxis.set_major_locator(plt.NullLocator())
    #print(f"x, y shape:: {np.array(x).shape} {y.shape}")
    if ys.ndim == 1:
        plt.plot(x, y, 'bo-')
        #plt.scatter(x, y)
    elif ys.ndim == 2:
        for i in range(len(Ylabels)):
            #plt.plot(x,ys[i,:], 'o-', label=Ylabels[i] )
            plt.plot(x,ys[i,:], label=Ylabels[i] )
    else:
        print(f"Error:: obscure in y-dim {ys.ndim}")
    ### ADD LEGEND
    plt.legend(loc=2)                   # locate after plot
    common_figure_after()              # comment to remove legend box 
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0

def draw_histogram(y, nbin, Lsave, fname):
    n, bins, patches = plt.hist(y, nbin)
    if Lsave:
        plt.savefig(fname, dpi=150)
    plt.show()

    return 0

def _mplot_2f(f1, f2, Lsave, figname, Title, Xtitle, Ytitle):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
    
def _f_draw(fname, dp, Lsave, figname, Title, Xtitle, Ytitle):
    x1=[]
    y1=[]
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    with open(fname, 'r') as f:
        for line in f:
            #if re.search("[a-zA-Z]", line):
            #    continue
            xy=line.split()
            if not xy[0]:
                del xy[0]
            xvalue=round(float(xy[0]), dp)
            yvalue=round(float(xy[1]), dp)

            x1.append(xvalue)
            y1.append(yvalue)
    print(x1, y1)
    #draw_2d(x1, y1, Lsave, fname )
    plt.plot(x1, y1)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0

def mplot_vector_one(x, y, Title=None, Xtitle=None, Ytitle=None ,Lsave=False):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
plot_line = mplot_vector_one

### this does not make an error in serial plotting
def mplot_vector_two(x, y, Title=None, Xtitle=None, Ytitle=None ,Lsave=False):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    for i in range(len(x)):
        print(x[i], y[i])
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0


def _mplot_2c(x, y, Lsave, figname, Title, Xtitle, Ytitle):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0

def _mplot_3c(x, y, y2, figname, Title, Xtitle, Ytitle, Lsave):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y )
    plt.plot(x, y2)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
def _mplot_2f3c(x, y1, y2,f1, f2, figname, Title, Xtitle, Ytitle, Lsave):
    #handles, labels = ax.get_legend_handels_labels()
    #ax.legend(handles, labels)
    plt.xlabel(Xtitle, fontsize=font_size)
    plt.ylabel(Ytitle, fontsize=font_size)
    plt.title(Title, fontsize=font_size)
    plt.plot(x, y1, label=f1)
    plt.plot(x, y2, label=f2)
    plt.legend()
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0


def barplot2(hash, name, title):
    """Makes a barplot of the fingerprint about the O atom."""
    fp = descriptor.fingerprints[hash][0]
    fig, ax = pyplot.subplots()
    ax.bar(range(len(fp[1])), fp[1])
    ax.set_title(title)
    ax.set_ylim(0., 2.)
    ax.set_xlabel('fingerprint')
    ax.set_ylabel('value')
    fig.savefig(name)

def barplot_y(ys, name=None, xlabel=None, ylabel=None, title=None):
    """Makes a barplot of the fingerprint about the O atom."""
    fig, ax = plt.subplots()
    ax.bar(range(1,len(ys)+1), ys)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    if name:
        fig.savefig(name)
