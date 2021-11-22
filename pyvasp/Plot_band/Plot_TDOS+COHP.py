#!/usr/bin/env python
import sys,time,glob,os
import matplotlib.pyplot as plt
import numpy as N
import pylab as plb
#import matplotlib.gridspec as gridspec



version = 20130814

nt = time.localtime()
now_time = "%s_%s_%s_%s%s%s" % (nt[0],nt[1],nt[2],nt[3],nt[4],nt[5])
usage = "Usage: %s [systemlabel] [systemlabel.EIG] [Pz.dat] [sigx.dat] [sigy.dat]" % sys.argv[0]
foottext = '\n Thank you\n## Kim,Hyo Seok (KAIST) <softmax1986@kaist.ac.kr>'

print "## Plotting PDOS+COHP"
print "## Version : %s \n" % version


def gen_pseudo_tickmark2(type_tick, loc, align, size=0.1):
    """
    Mimic tickmark just using a line...
    align = "On/Under the line"
    """
    import numpy as N
    pseudo_x = []; pseudo_y = []
    if type_tick == 'x':
        pseudo_y = N.arange(0, 2*size+0.005, 0.005)
        pseudo_x = N.zeros(len(pseudo_y)) + loc*N.ones(len(pseudo_y))
        if align == 'under':pseudo_y -= 2*size
    elif type_tick =='y':
        pseudo_x = N.arange(0, 2*size+0.01, 0.01)
        pseudo_y = N.zeros(len(pseudo_x)) + loc*N.ones(len(pseudo_x))
        if align == 'under':pseudo_x -= 2*size
    return pseudo_x, pseudo_y

#Check input file


if len(sys.argv) == 3:
    fermi_info = str(sys.argv[1]+'.EIG')
    T_dos = str(sys.argv[1]+'.DOS')
    cohp_1 = str(sys.argv[2]+'.cohp')

else :
    print usage
    print foottext
    sys.exit()
    

Pr_Gr = glob.glob('/home/softmax1986/Paper/Graphene_D_doped/7X7/Single_point_cal/Band_cal/Pr_Gr/Band_Pr_Gr.DOS')

f = open(Pr_Gr[0],"r")
f1 = open(fermi_info,"r")
f2 = open(T_dos,"r")
f3 = open(cohp_1, "r")

Fermi_Pr = -4.2668
Fermi = float(f1.readline())
Fermi_cohp = float(f3.readline())

Pr_energy = []
Pr_val = []

TDOS_energy = []
TDOS_up_val = []
TDOS_down_val = []
cohp_energy = []
cohp_d1_up = []
cohp_d1_down = []

TDOS_atoms = float(f2.readline())

for i in f.readlines():
    list_lines = i.split()
    Pr_energy.append(float(list_lines[0])-Fermi_Pr)
    Pr_val.append(float(list_lines[1])/98)

for i in f2.readlines():
    list_lines = i.split()
    TDOS_energy.append(float(list_lines[0])-Fermi)
    TDOS_up_val.append(float(list_lines[1])/TDOS_atoms)
    TDOS_down_val.append(float(list_lines[2])/TDOS_atoms)

for i in f3.readlines():
    list_lines = i.split()
    print list_lines
    cohp_energy.append(float(list_lines[0])-Fermi_cohp)
    cohp_d1_up.append(-float(list_lines[1])/124)        #number
    cohp_d1_down.append(-float(list_lines[2])/124)      #number

P_x = N.array(Pr_energy)
P_y = N.array(Pr_val)

T_x = N.array(TDOS_energy)
T_yu = N.array(TDOS_up_val)
T_yd = N.array(TDOS_down_val)

c_x = N.array(cohp_energy)
c1_yu = N.array(cohp_d1_up)
c1_yd = N.array(cohp_d1_down)

#plot data

fig = plt.figure(3, figsize=(15,15))
fig1 = fig.add_subplot(1.2,3,1,frame_on=1)
locs,labels = plb.yticks([],[])
locs,labels = plb.xticks([],[])
fig2 = fig.add_subplot(1.2,3,1.7,frame_on=1)
locs,labels = plb.yticks([],[])
locs,labels = plb.xticks([],[])
thickness = 3
frame  = fig1.get_frame()
frame2 = fig2.get_frame()
[i.set_linewidth(thickness) for i in fig1.spines.itervalues()]
[i.set_linewidth(thickness) for i in fig2.spines.itervalues()]
plb.rc('axes', linewidth=2)

#PDOS data Plot
fig1.plot(P_y,P_x,lw=3,color='grey',dashes=(15,10))
fig1.fill_betweenx(T_x,T_yd,0,color='orange')
fig1.fill_betweenx(T_x,T_yu,0,color='green',alpha=1)
fig1.axhline(0,color = 'black',lw=2,linestyle='--')

#COHP data Plot
fig2.plot(c1_yu, c_x, lw=5, color ='y',dashes=(5,5,5,5))
fig2.axhline(0,color = 'black',lw=2,linestyle='--')
fig2.axvline(0,color = 'black',lw=2)

# major pseudo xtickmarks #
for loc in [0.1]:
#    x, y = gen_pseudo_tickmark2('x', loc, 'under', size=0.03)
#    y += 2
#    figT1.plot(x,y,lw=1.5,color='black')
    x, y = gen_pseudo_tickmark2('x', loc, 'on', size=0.1)
    y -= 4
    fig1.plot(x,y,lw=3,color='black')

# major pseudo ytickmarks #
for loc in [-2,0,2]:
#    x, y = gen_pseudo_tickmark2('y', loc, 'under', size=0.09)
#    x -= 8
#    figT1.plot(x,y,lw=1.5,color='black')
    x, y = gen_pseudo_tickmark2('y', loc, 'on', size=0.001)
    x -=0 
    fig1.plot(x,y,lw=3,color='black')

for loc in [-0.02,0.02]:
#    x, y = gen_pseudo_tickmark2('x', loc, 'under', size=0.03)
#    y += 2
#    figT1.plot(x,y,lw=1.5,color='black')
    x, y = gen_pseudo_tickmark2('x', loc, 'on', size=0.1)
    y -= 4
    fig2.plot(x,y,lw=3,color='black')

fig1.axis([0,0.2,-4,4])
fig2.axis([-0.04,0.04,-4,4])
plt.savefig('PDOS+COHP.png',dpi=300)
plt.show()
