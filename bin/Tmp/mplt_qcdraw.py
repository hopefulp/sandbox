import matplotlib.pyplot as plt
### line draw function

XMIN=0
XMAX=3

YMIN=-0.5
YMAX=0.05           # for Ni-CO2 0.05, 1.0 for CO2 -charges
X_SHIFT=[0,1,2]
AB_GAP=0.4

Linewidth=2.0
Linewidth1=1.0

fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlim([XMIN,XMAX])
ax.set_ylim([YMIN,YMAX])

def f_draw(degeneracy, x, y, ehomo):
    if degeneracy == 1:
	    if float(y[0]) <= ehomo:
	        plt.plot(x, y, 'r' , lw=Linewidth)
	    else:
	        plt.plot(x, y, 'b', lw=Linewidth )
    elif degeneracy == 2:
        #print x, y
        x1=x[0:2]
        x2=x[2:4]
        if y[0] <= ehomo:
	        plt.plot(x1, y, 'r', x2, y, 'r', lw=Linewidth)
        else:
	        plt.plot(x1, y, 'b', x2, y, 'b', lw=Linewidth)
    else:
	    x1=x[0:2]
	    x2=x[2:4]
	    x3=x[4:6]
	    if y[0] <= ehomo:
	        plt.plot(x1, y, 'r', x2, y, 'r', x3, y, 'r', lw=Linewidth)
   	    else:
	        plt.plot(x1, y, 'b', x2, y, 'b', x3, y, 'b', lw=Linewidth)
    return 0

def f_draw_link(x, y):
    plt.plot(x, y, 'g', ls='dashed', lw=Linewidth1)

    return 0
### 
#key_word="Alpha MOs"
