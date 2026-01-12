"""
    2d plot library
    draw_2d
    fdraw   where f == "file"
    barplot
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import sys, re
import scipy
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
### choice between "import myplot_default"|call common_figurea

def lighten_color(color, ncolors):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    print(f"h l s {c[0]} {c[1]} {c[2]}")
    lmin = 0.2
    lmax = c[1]*1.9
    dl  = (lmax-lmin)/(ncolors-1)
    print(f"min, max, dl {lmin} {lmax} {dl}")
    color_rgb=[]
    for i in range(ncolors):
        print(f"lightness {lmin + dl*i}")
        color_rgb.append(colorsys.hls_to_rgb(c[0], lmin + dl*i, c[2]))
    amount=0.5 
    lightness = c[1]*0.3 + amount * (1 - c[1])
    print(f" lightness: {lightness }")
    #return colorsys.hls_to_rgb(c[0], lightness, c[2])
    return color_rgb
    #return colorsys.hls_to_rgb(c[0], amount * (1 - c[1]), c[2])

def lighten_2color(color, ncolors):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys

    color1=color[0]
    color2=color[1]
    nlight = int(ncolors/2)
    ndark  = ncolors -nlight

    c = colorsys.rgb_to_hls(*mc.to_rgb(color1))
    print(f"h l s {c[0]} {c[1]} {c[2]}")
    lmin = 0.25
    lmax = 0.75
    dl  = (lmax-lmin)/(ncolors-1)*2
    print(f"min, max, dl {lmin} {lmax} {dl}")
    color_rgb=[]
    color_light=[0.5, 0.8]
    for i in range(nlight):
        print(f"lightness {lmin + dl*i}")
        color_rgb.append(colorsys.hls_to_rgb(c[0], color_light[i], c[2]))
    c = colorsys.rgb_to_hls(*mc.to_rgb(color2))
    print(f"h l s {c[0]} {c[1]} {c[2]}")
    color_light=[0.8, 0.6, 0.3]
    for i in range(ndark):
        print(f"darkness {lmax - dl*i}")
        color_rgb.append(colorsys.hls_to_rgb(c[0], color_light[i], c[2]))
    #amount=0.5 
    #lightness = c[1]*0.3 + amount * (1 - c[1])
    #print(f" lightness: {lightness }")
    #return colorsys.hls_to_rgb(c[0], lightness, c[2])
    return color_rgb
    #return colorsys.hls_to_rgb(c[0], amount * (1 - c[1]), c[2])

def common_figure(ctype='dark', ncolor=4, Ltwinx=False):
    '''
    ctype   darken to change intensity
            cycle to use designated color turns
    '''
    if ctype == 'darken':
        import matplotlib.colors as ms
        import colorsys
    else:
        from cycler import cycler
    ### control figure size (2,6) for x-axis is 1/5
    fig = plt.figure(figsize=(10,6))         # def figsize=(10,6)
    ax = plt.axes()
    mpl.rcParams.update({'font.size':12})
    #ax.tick_params(axis='both', which='major', labelsize=25)
    #ax.tick_params(axis='x', labelsize=30)
    if ctype == 'darken':
        pass   
    else:
        if Ltwinx:
            print(f"make twinx axis")
            ax2 = ax.twinx()
            custom_cycler = (cycler(color=['r','m', 'orange']))
            custom_cycler2 = (cycler(color=['b', 'g']))
            ax.set_prop_cycle(custom_cycler)
            ax2.set_prop_cycle(custom_cycler2)
            return fig, ax, ax2
        else:
            if ncolor == 2:
                custom_cycler = (cycler(color=['r','g']))                                          # Figure 8(b)
                custom_cycler = (cycler(color=['darkcyan','b']))                                   # Figure S16(a)
                custom_cycler = (cycler(color=['r','darkcyan']))                                    # Figure S16(b)
            elif ncolor == 3:
                custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))       # Figure 8(a)
                custom_cycler = (cycler(color=['r','m','g','b']))
            elif ncolor == 4:
                custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))       # Figure 8(a)
                custom_cycler = (cycler(color=['r','m','g','b']))
            else:
                custom_cycler = (cycler(color=['r','m','g','b']))
            ax.set_prop_cycle(custom_cycler)
            return fig, ax
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


### twinx1: used for md.ene, normal data file
def mplot_twinx(x, y, iy_right, title=None, xlabel=None, ylabel=None, legend=None, Lsave=False, colors=None):
    '''
    called from "amp_plot_stat.py"
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    len(ylabel) == 2
    '''
    if iy_right:
        fig, ax, ax2 = common_figure(ctype='cycle', ncolor = len(y), Ltwinx=True)
    else:
        fig, ax = common_figure(ctype='cycle', ncolor = len(y))
    #fig, ax = plt.subplots(figsize=(12,8))
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
    ### try autocolor
    #plt.ylabel(ylabel1, fontsize=25, color='r')
    #ax.tick_params(axis='y', colors='r')
    ax.set_ylabel(ylabel1, fontsize=25)
    ax.tick_params(axis='y')
    ax.set_xlabel(xlabel, fontsize=25)
    #ax.set_ylim(-2,2)
    #ax.xaxis.set_major_locator(plt.NullLocator())
    #print(f"x, y shape:: {np.array(x).shape} {np.array(y).shape} and ylabel {ylabel} in {whereami()}")
    #ax2.set_ylabel(ylabel2, fontsize=25, color='g')
    ax2.set_ylabel(ylabel2, fontsize=25)
    #ax2.set_ylim(-2,2)
    tick_interval = 4.0
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    pls=[]
    print(f"y-right axis: index {iy_right} {len(ys)} {whereami()}")
    for i in range(len(ys)):
        if i in iy_right: 
            #if Colors:  color = Colors.pop(i)       #'tab:' + Colors.pop(0)
            #else:       color='tab:green'
            #plt.yticks(color='g')
            #p2, = ax2.plot(x, ys[i,:], '-', color='g', label=legend[i])
            print(f"{legend[i]} in y-right axis with color")
            if colors:
                p2, = ax2.plot(x, ys[i,:], '-',  label=legend[i], color=colors[i])
            else:
                p2, = ax2.plot(x, ys[i,:], '-',  label=legend[i])
            pls.append(p2)
        else:
            #ax2.tick_params(axis='y')
            #if Colors:  color = Colors.pop(i)       #color = 'tab:' + Colors.pop(0)
            #else:       color = 'tab:red'
            #plt.yticks(color='r')
            #p1, = ax.plot(x, ys[i,:], '-', color='r',  label=legend[i])
            if colors:
                p1, = ax.plot(x, ys[i,:], '-', label=legend[i], color=colors[i])
            else:
                p1, = ax.plot(x, ys[i,:], '-', label=legend[i])
            pls.append(p1)
    #plt.legend(pls, Ylabels, loc=2)
    ax.legend(loc=2)            # 2
    ax2.legend(loc=4)           # 1
    plt.legend()
    ### make the same y-scale in both size
    y_interval_ref = 18.30
    if iy_right:
        ymin, ymax = ax.get_ylim()
        y_interval = ymax - ymin
        if 'y_interval_ref' in locals():
            diff_interval = y_interval - y_interval_ref
            ymin_new = ymin + diff_interval/2.
            ymax_new = ymax - diff_interval/2.
            ax.set_ylim(ymin_new, ymax_new)
            y_interval = y_interval_ref
        ymin2, ymax2 = ax2.get_ylim()
        print(f"initial: ymin2 {ymin2} ymax2 {ymax2} y_interval {y_interval}")
        y2_interval = ymax2 - ymin2
        diff_interval = y2_interval - y_interval
        ymin2_new = ymin2 + diff_interval/2.
        ymax2_new = ymax2 - diff_interval/2.
        ax2.set_ylim(ymin2_new, ymax2_new)
        print(f"final: ymin2 {ymin2_new} ymax2 {ymax2_new}")
    #common_figure_after()
    #x_ticks = ['PP', 'PPP', 'PNP', 'PNP-bridged']
    #plt.xticks(x_ticks)
    #ax.set_xticklabels(x_ticks)
    #plt.locator_params(axis='x', nbins=10)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0

def make_double(x, y):
    xnew = []
    ynew = []
    line_pairs=[]
    for i, (x1, y1) in enumerate(zip(x, y)):
        xnew.append(x1)
        xnew.append(x1)
        ynew.append(y1)
        ynew.append(y1)
        if y1 is not None:
            line_pairs.append([[2*i, 2*i+1],[y1,y1]])
    return xnew, np.array(ynew), line_pairs


def mplot_levels(x, ys, title=None, xlabel=None, ylabel=None, legend=None,Lsave=False, colors=None):
    '''
    plot table 
    x   from 1st column
    y   other columns [0,1,2...] in  [[column],[size],...[size]]
    '''
    fig, ax = common_figure()
    ys = np.array(ys)
    if len(x) != ys.shape[1]:
        print(f"error in shape: x {len(x)}, ys.shape {ys.shape}")
    
    plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=35)
    plt.ylabel(ylabel, fontsize=35)
    
    if not colors:
        colors = ["blue", "red", "green", "orange", "purple"]
    else:
        colors = colors
    ### make levels: x=list, ys=array
    for i, y in enumerate(ys):
        xd, yd, linepair = make_double(x, y)
        print(f"Double: xd {xd} yd {yd}")
        xs = np.arange(len(xd))
        print(f"{np.max(xs)}")
        series = yd.astype(np.double)
        mask = np.isfinite(series)
        plt.plot(xs[mask], series[mask], '--', color=colors[i], label=legend[i])
        ### overdraw thick level
        for xl, yl in linepair: 
            plt.plot(xl,yl, '-', color=colors[i], lw=5 )
    #ax.yaxis.set_major_locator(ticker.MultipleLocator())
    #ax.set_xlim(-0.5, np.max(xs)+0.5)
    plt.legend(loc=1)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0


def auto_nvector(x,y):
    fig, ax = plt.subplots()
    for i in range(len(y)):
        plt.plot(x,y[i])
    plt.show()
    return 0

### For twinx, mplot_twinx
Lprint = 0
def mplot_nvector(x, y, plot_dict=None, Lsave=False, vertical=None, v_legend=None):
    '''
    input               python              numpy.shape
        x.shape         msize e.g. 18       (msize, )
        y.shape         (n, msize)          (n, msize)
    dx=1.0, input before
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    plot_dict   keys    xlabel, ylabel, xlim, title, colors
    '''
    print(f"input shape {np.array(x).shape} {np.array(y).shape}")

    ### get components of plot_dict
    if plot_dict:
        title   = plot_dict.get("title",    "")
        xlabel  = plot_dict.get("xlabel",   "")
        ylabel  = plot_dict.get("ylabel",   "")
        xlim    = plot_dict.get("xlim",     "")
        ylim    = plot_dict.get('ylim',     "")
        legends  = plot_dict.get("legend",   "")
        colors  = plot_dict.get("colors",   "")

    if 'colors' not in locals():
        fig, ax = common_figure(ncolor = len(legends))
        if Lprint: print("no color input")
    else:
        fig, ax = common_figure()
        if Lprint: print("colors")

    #plt.locator_params(axis='x', nbins=10)
    ys = np.array(y)
    ### x can be 1D or 2D
    xshape = np.array(x).shape
    if len(xshape) == 1:
        xdim = 1
        if len(x) != 0:
            xsize = len(x)
        else:
            xsize = ys.shape[0]
            x=range(xsize)
    else:
        xdim = 2
    print(f"xdim={xdim}")
    if xdim == 1:
        print(f"x= {xsize}, y={ys.shape} in {whereami()}()")
    else:
        print(f"x= {np.array(x).shape}, y={ys.shape} in {whereami()}()")
    plt.title(title)
    #xlabel = 'E - E$_F$ [eV]'
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    print(f"xlabel: {xlabel}, ys.ndim {ys.ndim} and loop for i = {len(ys)} ")
    '''
    print(f"legends {legends}")
    if not 'legends' in locals():
        legends=[]
        for i in range(len(ys)):
            legends.append(f'legend_{i:1d}')
    '''
    plt.rcParams['lines.linewidth'] = 2
    for i in range(len(ys)):
        ### treat xdim
        if xdim == 1:
            xs = x
        else:
            xs = x[i]
        ### is this for TDOS
        '''
        if re.match('T', legend[i]): # for what? for Tdos?
            d = scipy.zeros(len(ys[i,:]))
            print(f"shape d {np.array(d).shape}, legend ")
            #ax.fill_between(x, ys[i,:], where=ys[i,:]>=d, color=colors[i])
            plt.plot(xs,ys[i,:],  label=legend[i] , color=colors[i])
        else:
        '''
        #print(f"x dim {len(x)}, y dim {ys[i,:].shape} label {legend[i]}")
        print(f"x dim {len(x)}, y dim {ys[i,:].shape} in function {whereami()}()")
        #print(f"colors {colors}")
        if colors:
            ### if the final letter is 'f'
            if colors[i][-1] == 'f':
                color = colors[i][:-1]
                plt.fill_between(xs,ys[i,:], alpha=0.5, label=legends[i] , color=color)
            else:
                print(f"size x {len(xs)}, size y {len(list(ys[i,:]))}, size legend {len(legends)} size color {len(colors)}")
                plt.plot(xs,ys[i,:],  label=legends[i], color=colors[i])
        else:
            #if legends:
            if 0:
                plt.plot(xs,ys[i,:],  label=legends[i] )
            else:
                plt.plot(xs,ys[i,:])

    if vertical is not None:
        plt.axvline(x=vertical, color='k', linestyle='--', linewidth=1.5, label=v_legend)


    #else:
    #    print(f"Error:: obscure in y-dim {ys.ndim}")
    ### ADD LEGEND
    #plt.legend(loc=1)                # locate after plot
    if xlim:
        print(f"xlim in mplot2D.mplot_nvector {xlim}")
        plt.xlim(xlim)
    if ylim:
        #plt.ylim(bottom=0)
        plt.ylim(ylim)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    plt.show()
    if Lsave:
        fig.savefig("dos.png", dpi=150)
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
