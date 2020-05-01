"""
    2d plot library
    draw_2d
    fdraw   where f == "file"
    barplot
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
import sys
import my_chem 
import numpy as np
from myplot_default import *
#plt.switch_backend('agg')

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



def my_font(pack):
    if pack == 'amp':
        font={  'family': 'normal',
                'weight': 'bold',
                'size'  : 22 }
    else:
        print("package: %s is not included" % pack)
        sys.exit(0)
    mpl.rc('font', **font)
    return
### choice between "import myplot_default"|call common_figure
def common_figure():
    fig = plt.figure(figsize=(8,6))
    mpl.rcParams.update({'font.size':10})
    ax = plt.axes()
    ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(labelsize=25)

    custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))       # Figure 8(a)
    #custom_cycler = (cycler(color=['r','g']))                                          # Figure 8(b)
    #custom_cycler = (cycler(color=['darkcyan','b']))                                   # Figure S17(a)
    #custom_cycler = (cycler(color=['r','darkcyan']))                                    # Figure 8(b)
    #custom_cycler = (cycler(color=['orange','m','g','b']))
    ax.set_prop_cycle(custom_cycler)
    return fig, ax

def common_figure_after():
    plt.legend(loc=2)
    return 0


### draw_dots_two was upgraded to twinx
def draw_2subdots(y, h, title, suptitle, Ltwinx=None, escale=1.0,Colors=['r','b','o'], Ldiff=True):
    '''
    this makes error in serial plotting
    '''
    #fig, ax = common_figure()
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
def draw_dots_two(y, h, title, suptitle, Ltwinx=None, escale=1.0,Colors=['r','b','o'], Ldiff=True):
    '''
    this makes error in serial plotting
    '''
    #fig, ax = common_figure()
    escale = 1.0
    nlen = len(y)
    h_conv = np.array(h) * escale        # escale = my_chem.ev2kj
    y_conv = np.array(y) * escale
    ymin = min(y_conv)
    ymax = max(y_conv)
    y_width = ymax - ymin
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())
    max_res = abs(max(diff, key=abs))
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text_pos_x = nlen*0.85                  # 0.85, 0.2
    text_pos_y = max(y_conv)*0.2
    text="E_rms(test) = {:7.3f}\nE_maxres   = {:7.3f}".format(rmse, max_res)

    if Colors:  color = Colors.pop(0)       #'tab:' + Colors.pop(0)
    else:       color='tab:green'
    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    #mpl.rcParams.update({'font.size':22})
    plt.title(title, fontsize=size_title)
#    plt.suptitle(suptitle, x=0.5, y=0.92, va='top', fontsize=18)
    suptitle="\n"+suptitle
    plt.suptitle(suptitle, fontsize=size_tick)
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
        ax2.set_ylim(-0.5, 3.0)
        #ax2.plot(x, ys[1,:], 'o-', color=color)
        ax2.tick_params(axis='y', labelcolor='g', labelsize=size_tick)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    p1  = ax.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = ax.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    if Ltwinx:
        if Ldiff:
            p3, = ax2.plot(range(nlen), diff, c='g', label='difference')
            plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.45, 0.85))
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



def mplot_twinx(x, y, dx=1.0, Title=None, Xtitle=None, Ytitle=None, Ylabels=None, Lsave=False, Colors=None):
    '''
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[multi],[multi],...size]
    '''
    fig, ax = common_figure()
    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    #print(f"{x}-{ys}")

    plt.title(Title)
    if Xtitle:
        plt.xlabel(Xtitle, fontsize=15)
    plt.ylabel(Ytitle, fontsize=15)
    #ax.xaxis.set_major_locator(plt.NullLocator())
    #print(f"x, y shape:: {np.array(x).shape} {y.shape}")
    if Colors:  color = Colors.pop(0)       #color = 'tab:' + Colors.pop(0)
    else:       color = 'tab:red'
    ax.set_ylabel(Ylabels.pop(0), color=color)
    ax.plot(x, ys[0,:], 'o-', color=color)
    ax.tick_params(axis='y', labelcolor=color)

    ax2=ax.twinx()
    if Colors:  color = Colors.pop(0)       #'tab:' + Colors.pop(0)
    else:       color='tab:green'
    ax2.set_ylabel(Ylabels.pop(0), color=color)
    ax2.plot(x, ys[1,:], 'o-', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    #if ys.ndim == 1:
    #    plt.plot(x, y, 'bo-')
        #plt.scatter(x, y)
    #elif ys.ndim == 2:
    #    for i in range(len(Ylabels)):
    #        plt.plot(x,ys[i,:], 'o-', label=Ylabels[i] )
    #else:
    #    print(f"Error:: obscure in y-dim {ys.ndim}")

    #plt.legend(loc=2)                   # locate after plot
    #common_figure_after()              # legend box location
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0

def mplot_nvector(x, y, dx=1.0, Title=None, Xtitle=None, Ytitle=None, Ylabels=None, Lsave=False, Colors=None, Ltwinx=None):
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

def mplot_vector_one(y, Title=None, Xtitle=None, Ytitle=None ,Lsave=False):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
### this does not make an error in serial plotting
def mplot_vector_two(x, y, Title=None, Xtitle=None, Ytitle=None ,Lsave=False):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
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


def barplot(hash, name, title):
    """Makes a barplot of the fingerprint about the O atom."""
    fp = descriptor.fingerprints[hash][0]
    fig, ax = pyplot.subplots()
    ax.bar(range(len(fp[1])), fp[1])
    ax.set_title(title)
    ax.set_ylim(0., 2.)
    ax.set_xlabel('fingerprint')
    ax.set_ylabel('value')
    fig.savefig(name)

