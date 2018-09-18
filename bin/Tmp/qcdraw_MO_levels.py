#!/home/jackjack5/epd/bin/python
###  in developing for radical

import sys
import re
#from m_qcdraw import *
from mplt_qcdraw import *
from draw_MOs_mod import *
### read qchem output file
nfile=len(sys.argv)
if nfile <= 1 :
    print "Error with no arguments: input qchem outfile"
    exit(1)
qchem_outf=[]
for file in sys.argv:
    if re.search("out", file):
    	qchem_outf.append(file)
print "print 1: ", qchem_outf

### global variable 
#key_word="Alpha MOs"
key_word="Occupied"
key_lumo="Virtual"
key_beta="Beta MOs"
key_beta1="unrestricted"
end_word="-----"
fl_e_lumo=[]        # save lumo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
fl_line_ene=[]      # [ [eline1, eline2 ], [eline1, eline2...] [...
f_beta_tag=[0,0,0]  # [1, 0, 1] 1 for beta spin of  unrestricted calculation 
fid=0

Print_x="OFF"
for file in qchem_outf:
    e_lumo=100      # default for high energy
    tag_MO="NO"
    tag_lumo="NO"
    tag_beta="NO"
    inf=open(file, 'r')

    while 1:
        line=inf.readline()
        if not line: break 
 
        if tag_MO == "NO":
            if re.search(key_word, line):
                tag_MO="YES"
                line_ene=[]     # for job=optimization, whenever find new keyword, empty
                e_lumo_tmp=[]
        ### found MO
        else:
            #print line
            # if found virtual
            if re.search(key_lumo, line):
                tag_lumo="YES"
            # in case beta
            elif line == '\n':
                pass
            elif re.search(key_beta, line):
                tag_beta="YES" 
                tag_lumo="NO"   # for beta lumo
                f_beta_tag[fid]=1
                print "spin beta in ", file
            else:
                # if arrived the end line of MO energy
                if re.search(end_word, line):
                    tag_MO="NO"
                else:
                    # in case of beta, occupied appears again
                    if re.search(key_word, line):
                        pass
                    elif tag_lumo=="YES":
                        ene=re.split("\s+", line)
                        if ene[0]:
                            e_lumo=ene[0]
                        else:
                            e_lumo=ene[1]
                        tag_lumo="NO"
                        e_lumo_tmp.append(e_lumo)
                        line_ene.append(line)
                    else:
                    # if not anything, save
                        line_ene.append(line)
    # end of file
    inf.close()
    fid += 1
    if len(e_lumo_tmp) == 0:
	    e_lumo_tmp.append(0)
    if tag_beta == "YES":
	if len(e_lumo_tmp) == 1:
	    e_lumo_tmp.append(e_lumo_tmp[0])
    fl_e_lumo.append(e_lumo_tmp)
    fl_line_ene.append(line_ene)
### end of read file    
print "print 2: LUMO energy:: ", fl_e_lumo
fl_list_ene=[]      # [ [e1, e2... ], [e1, e2...] [...
for fline in fl_line_ene:
    list_ene=[]
    for eline in fline:
    	#print eline,
    	ene=re.split("\s+", eline)
    	for x in ene:
	    if x:       # remove empty entry
	    	list_ene.append(x)
    fl_list_ene.append(list_ene)
print "print 3: ",
for f_elist in fl_list_ene:
    print len(f_elist), f_elist


########## Draw mpl
### overall variable
fid=0
#print len(fl_list_ene)
for list_ene in fl_list_ene:
    ### initialize variables for each file drawing
    tag_draw_beta="NO"
    xlength=0.6
    x_3kinds=[]
    x_3kinds_b=[]
    ### for single level
    if f_beta_tag[fid] == 1:
        xlength /= 2
    x_l = XMIN + 0.2	+ X_SHIFT[fid]
    x_r = x_l  + xlength
    x_d = x_r - x_l	
    x_3kinds.append([x_l, x_r])
    if f_beta_tag[fid] == 1:
	    x_3kinds_b.append([x_l+AB_GAP, x_l+AB_GAP+x_d])
    ### for double degeneracy
    x_l1 = x_l		    	
    x_r1 = x_l + xlength/3  
    x_d  = x_r1-x_l1
    x_l2 = x_r1 + xlength/3 
    x_r2 = x_l2+x_d		  
    x2l=[x_l1, x_r1]
    x2r=[x_l2, x_r2]
    x2=x2l+x2r
    x_3kinds.append(x2)
    #print "print x2 interval: ", x2
    if Print_x=="ON":
        print "x2 range:: ", "%.2f "*len(x2) % tuple(x2)
    if f_beta_tag[fid] == 1:
        xbeta2l= [x_l1+AB_GAP, x_r1+AB_GAP]
        xbeta2r= [x_l2+AB_GAP, x_r2+AB_GAP]
        x2beta=xbeta2l+xbeta2r
        x_3kinds_b.append(x2beta)
        if Print_x=="ON":
            print "x2b range:: ", "%.2f "*len(x2beta) % tuple(x2beta)
    ### for triple degeneracy
    x_r1 = x_l + xlength/5
    x_d  = x_r1-x_l
    x_l2 = x_r1+0.05
    x_r2 = x_l2+x_d
    x_l3 = x_r2+0.05 
    x_r3 = x_l3+x_d
    x3l=[x_l1, x_r1]
    x3m=[x_l2, x_r2]
    x3r=[x_l3, x_r3]
    x3=x3l+x3m+x3r
    x_3kinds.append(x3)
    if Print_x=="ON":
        print "x3 range:: ", "%.2f "*len(x3) % tuple(x3)
    if f_beta_tag[fid] == 1:
        xbeta3l=[x_l1+AB_GAP, x_r1+AB_GAP]
        xbeta3m=[x_l2+AB_GAP, x_r2+AB_GAP]
        xbeta3r=[x_l3+AB_GAP, x_r3+AB_GAP]
        x3beta=xbeta3l+xbeta3m+xbeta3r
        x_3kinds_b.append(x3beta)
        if Print_x=="ON":
            print "x3b range:: ", "%.2f "*len(x3beta) % tuple(x3beta)
    i=0
    multi_degen=1
    e_lumo=float(fl_e_lumo[fid][0])
    #print "LUMO: ", e_lumo
    ### for H radical
    if len(list_ene) == 1:
	    y=[float(list_ene[0]), float(list_ene[0])]
	    f_draw(1, x_3kinds[0], y, e_lumo)
	    fid+=1
	    continue
    temp_ene=-1000
    tag_last_draw="NO"
    for v_ene in list_ene:
        if re.search("\*", v_ene): # in qchem outfile, 1st number is ******
            i+=1
            continue
        v_ene = float(v_ene)
    	if temp_ene==-1000:         # for initial energy value, just save 
            temp_ene=v_ene
            i+=1
            continue
    	y=[temp_ene, temp_ene]
    	##### if energies are same, increase multi degeneracy 
        if temp_ene == v_ene:
            multi_degen += 1
            ### if not last energy
            if i != len(list_ene)-1:
                i+=1
                continue
	        ### draw 3 plots if last ene and finish
            else:
                f_draw(multi_degen, x_3kinds[multi_degen-1], y, e_lumo)
		    #print "finish drawing"
	    ### if energies are different 
        else:
            ### draw 3 plots for previous energies
            f_draw(multi_degen, x_3kinds[multi_degen-1], y, e_lumo)
            # if last ene and no degenracy, plot additional energy
            if i == len(list_ene)-1:
                y=[v_ene, v_ene]
                f_draw(1, x_3kinds[0], y, e_lumo )
            # change into beta draw
            if(temp_ene > v_ene):
                print "Draw beta spin"
                tag_draw_beta="YES"
                e_lumo=float(fl_e_lumo[fid][1])
                x_3kinds=x_3kinds_b
        i+=1	
    	temp_ene=v_ene
    	multi_degen = 1
        ### end of one energy value
    fid+=1
    ### end of energy list
plt.show()


print "normal stop"	

