import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
### to use below: from myplot_default import *

### setting:: 2020.1.20 after system update, for AMP

#mpl.rcParams["figure.figsize"] = (8,6)          
mpl.rcParams.update({'font.size':15})           

#fig = plt.figure()
#ax = plt.axes()
#ax.tick_params(axis='both', which='major', labelsize=10)

### AMP in MLET
fig = plt.figure(figsize=(15,10))
ax = plt.axes()


custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))
ax.set_prop_cycle(custom_cycler)

size_title = 20
size_label = 18
size_tick = 15
text_x = 0.75
text_y = 0.8
text_twinx_x = 0.8
text_twinx_y = 0.9

#plt.title('AMP model error test-set')
#plt.ylabel('PE(kJ/mol)', fontsize=20)
plt.xlabel('data', fontsize=20)
#ax.tick_params(axis='both', which='major', labelsize=15)
#plt.scatter(x, y, 'r', x, y_bar, 'b')
#p1  = plt.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
#p2  = plt.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
#p3, = plt.plot(range(nlen), diff, label='difference')
#plt.plot(range(nlen), ones)
#plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'])



