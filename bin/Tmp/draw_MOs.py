#!/home/jackjack5/epd/bin/python
### 

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
end_word="-----"
i=0
fl_e_lumo=[]
fl_line_ene=[]
fid=0
for file in qchem_outf:
    e_lumo=100
    tag_MO="NO"
    tag_lumo="NO"
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
    	else:
	    if re.search(key_lumo, line):
	    	tag_lumo="YES"
	    else:
	    	if re.search(end_word, line):
		    tag_MO="NO"
	    	else:
		    if tag_lumo=="YES":
		    	ene=re.split("\s+", line)
		    	if ene[0]:
			    e_lumo=ene[0]
		    	else:
			    e_lumo=ene[1]
		    	tag_lumo="NO"
		    #print line,
		    line_ene.append(line)	    
    inf.close()
    fid += 1
    if e_lumo == 100:
	e_lumo=0
    fl_e_lumo.append(e_lumo)
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
ymin=-0.6
ymax=0.5
xlength=0.6
ax=plt.subplot(111)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin, ymax])
x_shift=[0,1,2]
j=0
#print len(fl_list_ene)
for list_ene in fl_list_ene:
    j+=1
    ### for single level
    x_l = xmin + 0.2 		+ x_shift[j-1]
    x_r = x_l  + xlength	
    x=[x_l, x_r]
    ### for double degeneracy
    x_l1 = x_l		    	
    x_r1 = x_l + xlength/3  
    x_l2 = x_r - xlength/3 
    x_r2 = x_r		  
    i=0
    multi_degen=1
    e_lumo=float(fl_e_lumo[j-1])
    print "LUMO: ", e_lumo
    ### for H radical
    if len(list_ene) == 1:
	y=[list_ene[0], list_ene[0]]
	plt.plot(x, y, 'r')
    for v_ene in list_ene:
	if re.search("\*", v_ene): ### in qchem outfile, 1st number is ******
	    continue
	v_ene = float(v_ene)
    	i += 1
    	if i == 1:
	    temp_ene=v_ene
	    continue
    	y=[temp_ene, temp_ene]
    	##### increase multi degeneracy 
    	if v_ene == temp_ene:
	    multi_degen += 1
	    if multi_degen >= 2:	## only for doubly degenerated
	    	#print "energy level has degeneracy of doublet"
	    	pass
	    else:
	    	print "error due to multi degeneracy"
    	##### if energy increases, draw
    	else:
	    if multi_degen == 1:
	    	if j==2:
		    print x, y, temp_ene, e_lumo
    	    	if temp_ene < e_lumo:
    	    	    plt.plot(x, y, 'r' )
    	    	else:
	    	    plt.plot(x, y, 'b' )
	    elif multi_degen >= 2:
	    	x2l=[x_l1, x_r1]
	    	x2r=[x_l2, x_r2]
		#if j == 2:
		#    print temp_ene, e_lumo
	    	if temp_ene < e_lumo:
	    	    plt.plot(x2l, y, 'r', x2r, y, 'r')
	    	else:
		    plt.plot(x2l, y, 'b', x2r, y, 'b')
		if multi_degen > 2:
		    print "Be cautious:: triply degenerated"
	    	multi_degen = 1
	    else:
	    	print "Error in drawing triple level"
	    	exit(10)
    	    temp_ene=v_ene
	### for H2 with 1 lumo
	if i == len(list_ene):
	    y=[temp_ene, temp_ene]
	    plt.plot(x, y, 'b' )   # for its virtual
plt.show()


print "normal stop"	

