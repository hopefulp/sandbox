import matplotlib as mpl
import matplotlib.pyplot as plt
#from mplt_fig_set import *
#from qcout_ini import *
import os 
#### print option
tag_print=1

### line draw function
TEMP_ENE=-1000
print_draw=0

#### MatPlotLib ini #####
#plt.style.use('classic')
mpl.rcParams['xtick.labelsize']=16
mpl.rcParams['ytick.labelsize']=16
mpl.rcParams['lines.linewidth']=5
mpl.rcParams['figure.figsize']=(15, 12)
mpl.rcParams['axes.titlesize']= 'large'
#plt.ylabel('E(Hr)')
## arrow


### mplot_qcdraw 
### x divides 3 parts for 3 files
XMIN=0
XMAX=1
XLength=0.07
XL_2=XLength/2.
# for Ni-CO2 0.05, 1.0 for CO2 -charges
YMAX    =  0.61
YMIN    =  -0.55

"""
if re.search("1-P", files[0]):
    YMAX    =  0.05
    YMIN    = -0.4
elif re.search("2-P", files[0]):
    YMAX    =  0.05
    YMIN    = -0.3
elif re.search("3-P", files[0]):   
    YMAX    =  0.2
    YMIN    = -0.3
else:
    YMAX    = 0.5       # 0.2
    YMIN    = -0.6      # -0.7
"""


#### main module cannot give a value to sub-module
#### another module is required to get Y values

X_SHIFT = [0,1,2,3,4]
AB_GAP  = 0.4

Linewidth=2.0
Linewidth1=1.0

fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlim([XMIN,XMAX])
ax.set_ylim([YMIN,YMAX])

def mplot_arrow_homo(x, y, l_split_s, spin):
    """
        x, y is 2-element rank=1
    """        
    if l_split_s==0:
        xm1=(x[0]+x[1])/2.-0.02
        xm2=(x[0]+x[1])/2.+0.02
        plt.arrow(xm1, y[0]-0.02, 0.0, 0.04, width=0.002, length_includes_head=True)
        plt.arrow(xm2, y[0]+0.02, 0.0, -0.04, width=0.002, length_includes_head=True)
    else:
        xm1=(x[0]+x[1])/2.
        if spin==0:
            plt.arrow(xm1, y[0]-0.02, 0.0, 0.04, width=0.002, length_includes_head=True)
        else:
            plt.arrow(xm1, y[0]+0.02, 0.0, -0.04, width=0.002, length_includes_head=True)
        
    return 0

def mplt_level(degeneracy, x, y, ehomo, s_tag, abspin):
    """
    mplt_level arges::  degeneracy, x value list, y value list, homo energy
    draw one level with a few lines depending on degeneracy
    degeneracy depends on the size of x_list
    maximum degeneracy is 3
    """
    if degeneracy == 1:
        if float(y[0]) <= ehomo:
            #print("degeneracy-1", x[0], y[0])
            plt.plot(x, y, 'r' , lw=Linewidth)
            ### in case: only homo is drawn by arrow
            #if y[0] == ehomo:
            mplot_arrow_homo(x,y,s_tag,abspin)
        else:
            plt.plot(x, y, 'b', lw=Linewidth )
    elif degeneracy == 2:
        print (x, y)
        x1=x[0:2]
        x2=x[2:4]
        if y[0] <= ehomo:
            plt.plot(x1, y, 'r', x2, y, 'r', lw=Linewidth)
            #if y[0] == ehomo:
            mplot_arrow_homo(x1,y,s_tag,abspin)
            mplot_arrow_homo(x2,y,s_tag,abspin)
        else:
            plt.plot(x1, y, 'b', x2, y, 'b', lw=Linewidth)
    elif degeneracy == 3:
        x1=x[0:2]
        x2=x[2:4]
        x3=x[4:6]
        if y[0] <= ehomo:
            plt.plot(x1, y, 'r', x2, y, 'r', x3, y, 'r', lw=Linewidth)
            mplot_arrow_homo(x1,y,s_tag,abspin)
            mplot_arrow_homo(x2,y,s_tag,abspin)
            mplot_arrow_homo(x3,y,s_tag,abspin)
        else:
            plt.plot(x1, y, 'b', x2, y, 'b', x3, y, 'b', lw=Linewidth)
    else:
        print ("Error in n_degeneracy in function mplt_level")
        exit(55)
    return 0

def mplt_level_link(x, y):
    plt.plot(x, y, 'g', ls='dashed', lw=Linewidth1)
    #print x, y, " in ", os.path.basename(__file__)
    return 0
### 
#key_word="Alpha MOs"

#### xmin, xmax for a file
def fx_region(x0, x1, n_f, i_f):
    """
        returns x_min, x_max - X diveded by number of files
    """
    dx = float(x1-x0)/n_f
    xmin = i_f * dx
    xmax = (i_f+1) * dx
    return xmin, xmax

def Xrange_nf_fixed_x_length(nfile,ifile,L_beta,ab, degen):
    """
    Xrange_nf args:: n_files, index_file, only_alpha(0)_or_beta_exists(1), alpha(0)_or_beta(1), degeneracy
    xrange: 0 ~ to XMAX
    nfile: divide xrange with number of files
    ifile: indicates each section of xrange by xmin, xmax of each file
    beta: xmin~xmax is divided by half
    alpha locates 1st half and beta locates 2nd half
        dx: half of each file is divided by degeneracy
    """
    #### xrange for files
    xrange=float(XMAX-XMIN)/nfile
    #### devide xrange into partsf for degeneracy
    #npart=2*degen+1
    #dx=float(xrange)/float(npart)   # devide Xrange into parts

    #### xmin, xmax for a file
    if_xmin, if_xmax = fx_region(XMIN, XMAX, nfile, ifile)
    #print "xrange", xmin, xmax
    ### xmin, xmax for alpha or beta in the file            
    if L_beta:
        x_2=float(if_xmax-if_xmin)/2
        if ab==0:
            if_xmax=if_xmin+x_2
        elif ab==1:
            if_xmin+=x_2
        else:
            print ("error in function draw")
            exit(97)
        #dx=float(x_2)/float(npart)
    if_x_center = (if_xmax - if_xmin) /2. + if_xmin

    if degen == 1:
        x1 = if_x_center - XL_2
        x2 = if_x_center + XL_2
        return x1, x2
    elif degen == 2:
        x1 = if_x_center - XL_2*3
        x2 = if_x_center - XL_2*1
        x3 = if_x_center + XL_2*1
        x4 = if_x_center + XL_2*3
        return x1, x2, x3, x4
    elif degen==3:
        x1 = if_x_center - XL_2*4
        x2 = if_x_center - XL_2*2
        x3 = if_x_center - XL_2*1
        x4 = if_x_center + XL_2*1
        x5 = if_x_center + XL_2*2
        x6 = if_x_center + XL_2*4
        return x1, x2, x3, x4, x5, x6
    elif degen==4:
        x1 = if_x_center - XL_2*5
        x2 = if_x_center - XL_2*3
        x3 = if_x_center - XL_2*1
        x4 = if_x_center + XL_2*1
        x5 = if_x_center + XL_2*3
        x6 = if_x_center + XL_2*5
        return x1, x2, x3, x4, x5, x6, x7, x8
        
    return 0

def Xrange_nf_cal_dx(nfile,ifile,L_beta,ab, degen):
    """
    Xrange_nf args:: n_files, index_file, only_alpha(0)_or_beta_exists(1), alpha(0)_or_beta(1), degeneracy
    xrange: 0 ~ to XMAX
    nfile: divide xrange with number of files
    ifile: indicates each section of xrange by xmin, xmax of each file
    beta: xmin~xmax is divided by half
    alpha locates 1st half and beta locates 2nd half
        dx: half of each file is divided by degeneracy
    """
    #### xrange for files
    xrange=float(XMAX-XMIN)/nfile
    #### devide xrange into partsf for degeneracy
    npart=2*degen+1
    dx=float(xrange)/float(npart)   # devide Xrange into parts

    #### xmin, xmax for a file
    xmin, xmax = fx_region(XMIN, XMAX, nfile, ifile)
    #print "xrange", xmin, xmax
    ### xmin, xmax for alpha or beta in the file            
    if L_beta:
        x_2=float(xmax-xmin)/2
        if ab==0:
            xmax=xmin+x_2
        elif ab==1:
            xmin+=x_2
        else:
            print ("error in function draw")
            exit(97)
        dx=float(x_2)/float(npart)

    x1=xmin+dx
    x2=x1+dx
    x3=x2+dx
    x4=x3+dx
    x5=x4+dx
    x6=x5+dx
    x7=x6+dx
    x8=x7+dx
    if degen == 1:
        return x1, x2
    elif degen == 2:
        return x1, x2, x3, x4
    elif degen==3:
        return x1, x2, x3, x4, x5, x6
    elif degen==4:
        return x1, x2, x3, x4, x5, x6, x7, x8
        
    return 0

def mplt_ab_draw(nfile, fid, beta_tag, ab_tag, diction, e_homo):
    """
    type ab_draw(int, int, int, int, hash(int, float), float)
    args:: n_files, index_file, beta_exists?, alpha_or_beta, imo_ene_dictionary, homo_energy
        sub function: Xrange_nf, mplt_level
        degeneracy is calculated using dictionary
    """
    #### change into sorted dictionary: but list of (key, value) tuples
    if tag_print >= 2:
        print ("func: ab_draw ", e_homo)
    l_dic=sorted(diction.items())     
    nkeys=len(l_dic)

    #### initialize variables for each file drawing
    #tag_draw_beta="NO"
    #### MO level ragne is classified by nfiles
    i=0
    degen_eline=1
    #e_homo=float(FL_homo_ene[fid][0])
    #### for H radical
    if nkeys == 1:
        y=[l_dic[0][1], l_dic[0][1]]
        mplt_level(1, x_3kinds[0], y, e_homo,0, 0)
        fid+=1
        #continue
    temp_ene=TEMP_ENE
    tag_last_draw="NO"

    #if re.search("\*", v_ene): # in qchem outfile, 1st low number might be ******
    #### When use dictionary, be needed to sort keys 
    # fix the size of bar
    for tup in l_dic:
        i+=1
        v_ene=tup[1]
        if print_draw:
            print (v_ene)
        if temp_ene==TEMP_ENE:         # for initial energy value, just save 
            temp_ene=v_ene
            continue
        #### draw the previous energy y
        y=[temp_ene, temp_ene]
        #### if energies are same, increase multi degeneracy 
        if temp_ene == v_ene:
            degen_eline += 1
            ### if not last energy
            if i != nkeys:
                continue
            ### draw 3 plots if last ene and finish
            else:
                x_r=Xrange_nf_fixed_x_length(nfile, fid, beta_tag, ab_tag, degen_eline)
                mplt_level(degen_eline, x_r, y, e_homo, beta_tag, ab_tag)
            #print "finish drawing"
        ### Draw if energies are different 
        else:
            ### draw 3 kind of plots with degeneracy for previous energies
            if degen_eline > 4:
                print ("ERROR:: so many degeneracy in drawing")
                exit(10)
            #print "degen_eline = ", degen_eline
            #print nfile, fid, beta_tag, ab_tag
            x_r=Xrange_nf_fixed_x_length(nfile, fid, beta_tag, ab_tag, degen_eline)
            #print "ab_tag: ", ab_tag,  x_r
            mplt_level(degen_eline, x_r, y, e_homo , beta_tag, ab_tag)
            # if last ene and no degenracy, plot additional energy
            if i == nkeys:
                y=[v_ene, v_ene]
                #print nfile, fid, beta_tag, ab_tag
                #print Xrange_nf(nfile, fid, beta_tag, ab_tag, 1)
                x_r=Xrange_nf_fixed_x_length(nfile, fid, beta_tag, ab_tag, 1)
                #print x_r
                mplt_level(1, x_r, y, e_homo, beta_tag, ab_tag )
        temp_ene=v_ene
        degen_eline = 1
        #### end of one energy value
    #### end of energy list
    return 0


def main():
    print("matplotlib.rcParams are set")
    print("def mplt_ab_draw:: main function for draw including alpha, beta")
    print("            calls: Xrange_nf_fixed_x_length for x dividing")
    print("def Xrange_nf_fixed_x_length:: x's for 5 files")


if __name__ == "__main__":
    main()
