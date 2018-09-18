#!/home/jackjack5/epd/bin/python
###  in developing for radical

import sys
import re

nfile=len(sys.argv)
if nfile <= 1 :
    print "Error with no arguments: input qchem outfile"
    exit(1)
qchem_outf=[]
for file in sys.argv:
    if re.search("out", file):
    	qchem_outf.append(file)

print qchem_outf
#key_word="Alpha MOs"
key_word="Occupied"
key_lumo="Virtual"
key_beta="Beta MOs"
key_beta1="unrestricted"
end_word="-----"
i=0
fl_e_lumo=[]
fl_line_ene=[]
f_beta_tag=[0,0,0]
fid=0
for file in qchem_outf:
    e_lumo=100
    tag_MO="NO"
    tag_lumo="NO"
    tag_beta="NO"
    inf=open(file, 'r')

    while 1:
    	line=inf.readline()
    	if not line: break 
 
    	if tag_MO == "NO":
    	    if re.search(key_word, line):
	    	#print line,
	    	#print "have found"
	    	tag_MO="YES"
	    	line_ene=[]
    		e_lumo_tmp=[]
    	else: # after Occupied
	    #print line
	    if re.search(key_lumo, line):
	    	tag_lumo="YES"
	    elif line == '\n':
		pass
	    elif re.search(key_beta, line):
		tag_beta="YES" 
		tag_lumo="NO"
		f_beta_tag[fid]=1
		print "here is beta"
	    else:
	    	if re.search(end_word, line):
		    tag_MO="NO"
	    	else:
		    if re.search(key_word, line):
			pass
		    # to save e_lumo
		    elif tag_lumo=="YES":
		    	ene=re.split("\s+", line)
		    	if ene[0]:
			    e_lumo=ene[0]
		    	else:
			    e_lumo=ene[1]
		    	tag_lumo="NO"
		    	line_ene.append(line)
			e_lumo_tmp.append(e_lumo)
		    else:
			line_ene.append(line)
		    #print line,
    inf.close()
    fid += 1
    if len(e_lumo_tmp) == 0:
	e_lumo_tmp.append(0)
    if tag_beta == "YES":
	if len(e_lumo_tmp) == 1:
	    e_lumo_tmp.append(e_lumo_tmp[0])
    fl_e_lumo.append(e_lumo_tmp)
    fl_line_ene.append(line_ene)
print "LUMO energy:: ", fl_e_lumo
fl_list_ene=[]
for fline in fl_line_ene:
    list_ene=[]
    for eline in fline:
    	print eline,
    	ene=re.split("\s+", eline)
    	for x in ene:
	    if x:
	    	list_ene.append(x)
    fl_list_ene.append(list_ene)
print fl_list_ene
########## Draw mpl
import matplotlib.pyplot as plt

fig=plt.figure()
#ax=plt.axes(xlim(0,2))
xmin=0
xmax=3
ymin=-1.6
ymax=0.5
xlength=0.6
ax=plt.subplot(111)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin, ymax])
x_shift=[0,1,2]
fid=0
tag_draw_beta="NO"
ab_gap=0.4
#print len(fl_list_ene)
for list_ene in fl_list_ene:
    ### for single level
    if f_beta_tag[fid] == 1:
	xlength /= 2
    x_l = xmin + 0.2	+ x_shift[fid]
    x_r = x_l  + xlength
    x_d = x_r - x_l	
    x=[x_l, x_r]
    if f_beta_tag[fid] == 1:
	xbeta=[x_l+ab_gap, x_l+ab_gap+x_d]	
    ### for double degeneracy
    x_l1 = x_l		    	
    x_r1 = x_l + xlength/3  
    x_d  = x_r1-x_l1
    x_l2 = x_r1 + xlength/3 
    x_r2 = x_l2+x_d		  
    x2l=[x_l1, x_r1]
    x2r=[x_l2, x_r2]
    if f_beta_tag[fid] == 1:
	xbeta2l= [x_l1+ab_gap, x_r1+ab_gap]
	xbeta2r= [x_l2+ab_gap, x_r2+ab_gap]
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
    print x3l, x3m, x3r
    if f_beta_tag[fid] == 1:
	xbeta3l=[x_l1+ab_gap, x_r1+ab_gap]
	xbeta3m=[x_l2+ab_gap, x_r2+ab_gap]
	xbeta3r=[x_l3+ab_gap, x_r3+ab_gap]
	print xbeta3l, xbeta3m, xbeta3r

    i=0
    multi_degen=1
    e_lumo=float(fl_e_lumo[fid][0])
    print "LUMO: ", e_lumo
    ### for H radical
    if len(list_ene) == 1:
	y=[list_ene[0], list_ene[0]]
	plt.plot(x, y, 'r')
    temp_ene=-1000
    tag_last_draw="NO"
    for v_ene in list_ene:
	if re.search("\*", v_ene): ### in qchem outfile, 1st number is ******
	    i+=1
	    continue
	v_ene = float(v_ene)
    	if temp_ene==-1000:
	    temp_ene=v_ene
	    i+=1
	    continue
    	y=[temp_ene, temp_ene]
    	##### increase multi degeneracy 
    	if temp_ene == v_ene and i != len(list_ene)-1:
	    #print i, len(list_ene), multi_degen, v_ene
	    multi_degen += 1
	    i+=1
	    continue
	else:
	    # if last ene and degeneracy
	    if temp_ene == v_ene:
		multi_degen += 1
    	    ##### if energy increases, draw
	    print "draw previous energy", i, v_ene
	    if multi_degen == 1:
    	    	if temp_ene < e_lumo:
    	    	    plt.plot(x, y, 'r' )
    	    	else:
	    	    plt.plot(x, y, 'b' )
	    elif multi_degen == 2:
	    	if temp_ene < e_lumo:
	    	    plt.plot(x2l, y, 'r', x2r, y, 'r')
	    	else:
		    plt.plot(x2l, y, 'b', x2r, y, 'b')
		if i == len(list_ene)-1:
		    tag_last_draw="YES"
	    elif multi_degen == 3:
		if tag_draw_beta == "YES":
		    print "why here? draw beta"
		else:
		    print "draw alpha"
	    	if temp_ene < e_lumo:
	    	    plt.plot(x3l, y, 'r', x3m, y, 'r', x3r, y, 'r')
	    	else:
		    plt.plot(x3l, y, 'b', x3m, y, 'b', x3r, y, 'b')
		    print x3l, x3m, x3r
		if i == len(list_ene)-1:
		    tag_last_draw="YES"
	    # if last ene and no degeneracy 
	    if i == len(list_ene)-1 and tag_last_draw=="NO":
	    	y=[v_ene, v_ene]
	    	plt.plot(x, y, 'b' )
	    if(temp_ene > v_ene):
		print "x vector changed into beta"
		tag_draw_beta="YES"
		x=xbeta
		x2l=xbeta2l
		x2r=xbeta2r
		x3l=xbeta3l
		x3m=xbeta3m
		x3r=xbeta3r
	i+=1	
    	temp_ene=v_ene
    	multi_degen = 1
        ### end of one energy value
    fid+=1
    ### end of energy list
print "multi = ", multi_degen
plt.show()
print "normal stop"	

