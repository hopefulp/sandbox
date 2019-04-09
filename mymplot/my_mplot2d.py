"""
    2d plot library
    function
        draw_2d
        fdraw   where f == "file"
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import my_chem 
import numpy as np

### do not open figure here
#fig = plt.figure(figsize=(15,10))
#ax = plt.axes()

#ax = plt.subplot(111)
#ax = fig.gca()

### font size
#font = {'size': 15}

#mpl.rcParams['legend.fontsize'] = 10
#mpl.rc('font', **font)
#font_size=15
#mpl.rcParams.update({'font.size':22})

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

def common_figure():
    fig = plt.figure(figsize=(15,10))
    ax = plt.axes()
    mpl.rcParams.update({'font.size':22})
    ax.tick_params(axis='both', which='major', labelsize=15)
    return fig, ax

def draw_dots_two(y, h, title, suptitle):
    '''
    this makes error in serial plotting
    '''
    fig = plt.figure(figsize=(15,10))
    ax = plt.axes()

    nlen = len(y)
    h_conv = np.array(h) * my_chem.ev2kj
    y_conv = np.array(y) * my_chem.ev2kj
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())
    max_res = abs(max(diff, key=abs))
    #max_res = max(diff, key=abs)
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text_pos_x = nlen*0.85
    text_pos_y = max(y_conv)*0.2
    text="E_rms(test) = {:7.3f}\nE_maxres = {:7.3f}".format(rmse, max_res)
    plt.text(text_pos_x, text_pos_y,text, fontsize=20)


    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    mpl.rcParams.update({'font.size':22})
    plt.title(title, fontsize=40)
#    plt.suptitle(suptitle, x=0.5, y=0.92, va='top', fontsize=18)
    plt.suptitle(suptitle, fontsize=18)
    plt.ylabel('PE(kJ/mol)', fontsize=20)
    plt.xlabel('data', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    p1  = plt.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = plt.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    p3, = plt.plot(range(nlen), diff, label='difference')
    plt.plot(range(nlen), ones)
    plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
    plt.show()
    return rmse
def xtitle_font(tit):
    st = "\'{}\', fontsize=20".format(tit)
    print(st)
    return st

def mplot_nvector(x, y, Title=None, Xtitle=None, Ytitle=None, Ylabels=None, Lsave=False):
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
    #print("hmm: {}".format(xsize))
    plt.xticks(np.arange(min(x), max(x)+1, int(max(x)/10)))

    plt.title(Title)
    plt.xlabel(Xtitle, fontsize=20)
    plt.ylabel(Ytitle, fontsize=20)
    if ys.ndim == 1:
        #plt.plot(x, y)
        plt.scatter(x, y)
    elif ys.ndim == 2:
        for i in range(len(Ylabels)):
            plt.plot(x,ys[:,i], label=Ylabels[i])

    plt.legend(loc=1)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)


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
 
