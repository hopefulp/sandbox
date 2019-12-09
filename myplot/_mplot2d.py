"""
    2d plot library
    function
        draw_2d
        fdraw   where f == "file"
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

### font size
font = {'size': 15}

mpl.rcParams['legend.fontsize'] = 10
mpl.rc('font', **font)
font_size=15

fig = plt.figure()
#ax = plt.subplot(111)
ax = fig.gca()

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
#def draw_2d(x, y, Lsave, fname):
#    plt.plot(x, y)
#    plt.show()
#    return 0
    
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
 
