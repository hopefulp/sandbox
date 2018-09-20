"""
    2d plot library
    function
        draw_2d
        fdraw   where f == "file"
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['legend.fontsize'] = 20
fig = plt.figure()
ax = plt.subplot(111)

def draw_histogram(y, nbin, Lsave, fname):
    n, bins, patches = plt.hist(y, nbin)
    if Lsave:
        plt.savefig(fname, dpi=150)
    plt.show()

    return 0

#def draw_2d(x, y, Lsave, fname):
#    plt.plot(x, y)
#    plt.show()
#    return 0
    
def f_draw(fname, dp, Lsave, figname, Title, Xtitle, Ytitle):
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

