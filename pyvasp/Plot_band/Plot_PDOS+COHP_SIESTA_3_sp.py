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


if len(sys.argv) == 8:
    fermi_info = str(sys.argv[1]+'.EIG')
    T_dos = str(sys.argv[2]+'.dat')
    sigx_dos = str(sys.argv[3]+'.dat')
    sigy_dos = str(sys.argv[4]+'.dat')
    cohp_1 = str(sys.argv[5]+'.cohp')
    cohp_2 = str(sys.argv[6]+'.cohp')
    cohp_3 = str(sys.argv[7]+'.cohp')

else :
    print usage
    print foottext
    sys.exit()
    

Pr_Gr = glob.glob('/home/softmax1986/Paper/Graphene_D_doped/7X7/Band_cal/Pr_Gr/Pr_Pz.dat')

f = open(Pr_Gr[0],"r")
f1 = open(fermi_info,"r")
f2 = open(T_dos,"r")
f3 = open(sigx_dos,"r")
f4 = open(sigy_dos,"r")
f5 = open(cohp_1, "r")
f6 = open(cohp_2, "r")
f7 = open(cohp_3, "r")

Fermi_Pr = -4.2670
Fermi = float(f1.readline())
#Fermi_cohp = float(f4.readline())

Pr_energy = []
Pr_val = []

TDOS_energy = []
TDOS_up_val = []
TDOS_down_val = []
sigx_up_val = []
sigx_down_val = []
sigy_up_val = []
sigy_down_val = []
cohp_energy = []
cohp_d1_up = []
cohp_d1_down = []
cohp_d2_up = []
cohp_d2_down = []
cohp_d3_up = []
cohp_d3_down = []

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
    sigx_up_val.append(float(list_lines[1])/TDOS_atoms)
    sigx_down_val.append(float(list_lines[2])/TDOS_atoms)

for i in f4.readlines():
    list_lines = i.split()
    sigy_up_val.append(float(list_lines[1])/TDOS_atoms)
    sigy_down_val.append(float(list_lines[2])/TDOS_atoms)

for i in f5.readlines():
    list_lines = i.split()
    cohp_energy.append(float(list_lines[0]))
    cohp_d1_up.append(float(list_lines[1]))
    cohp_d1_down.append(float(list_lines[2]))

for i in f6.readlines():
    list_lines = i.split()
    cohp_d2_up.append(float(list_lines[1]))
    cohp_d2_down.append(float(list_lines[2]))

for i in f7.readlines():
    list_lines = i.split()
    cohp_d3_up.append(float(list_lines[1]))
    cohp_d3_down.append(float(list_lines[2]))

P_x = N.array(Pr_energy)
P_y = N.array(Pr_val)

T_x = N.array(TDOS_energy)
T_yu = N.array(TDOS_up_val)
T_yd = N.array(TDOS_down_val)

sig_yu = ((N.array(sigx_up_val)+N.array(sigy_up_val))/2)
sig_yd = ((N.array(sigx_down_val)+N.array(sigy_down_val))/2)

c_x = N.array(cohp_energy)
c1_yu = N.array(cohp_d1_up)
c1_yd = N.array(cohp_d1_down)
c2_yu = N.array(cohp_d2_up)
c2_yd = N.array(cohp_d2_down)
c3_yu = N.array(cohp_d3_up)
c3_yd = N.array(cohp_d3_down)

#plot data

fig = plt.figure(3, figsize=(18,18))
fig1 = fig.add_subplot(1,2,1,frame_on=1)
locs,labels = plb.yticks([],[])
locs,labels = plb.xticks([],[])
fig2 = fig.add_subplot(1,4,2.43,frame_on=1)
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
fig1.fill_betweenx(T_x,T_yu,0,color='red',alpha=0.8)
fig1.fill_betweenx(T_x,sig_yd,0,color='chartreuse',alpha=1)
fig1.fill_betweenx(T_x,sig_yu,0,color='g',alpha=1)
fig1.plot(T_yu,T_x,lw=2,color='black')
fig1.plot(sig_yu,T_x,lw=2,color='black')
fig1.plot(T_yd,T_x,lw=2,color='black',dashes=(15,10,5,10))
fig1.plot(sig_yd,T_x,lw=2,color='black',dashes=(15,10,5,10))
fig1.axhline(0,color = 'black',lw=3,linestyle='--')

#COHP data Plot
fig2.plot(c1_yu, c_x, lw=6, color ='blue')
fig2.plot(c1_yd, c_x, lw=6, color ='blue',dashes=(15,10))
fig2.plot(c2_yu, c_x, lw=6, color ='red')
fig2.plot(c2_yd, c_x, lw=6, color ='red',dashes=(15,10))
fig2.plot(c3_yu, c_x, lw=6, color ='chartreuse')
fig2.plot(c3_yd, c_x, lw=6, color ='chartreuse',dashes=(15,10))
fig2.axhline(0,color = 'black',lw=3,linestyle='--')
fig2.axvline(0,color = 'black',lw=3)

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

fig1.axis([0,0.15,-4,4])
fig2.axis([-0.1,0.1,-4,4])
plt.savefig('PDOS+COHP.png',dpi=300)
plt.show()
