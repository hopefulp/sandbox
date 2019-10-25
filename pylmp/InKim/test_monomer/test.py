#!/home/noische/python

from monomer import *

g=RandomPolymer()
g.Dendron="/qcfs/noische/research/PEI/201510/tertiary.b3lyp.final.bgf"
g.Linear="/qcfs/noische/research/PEI/201510/secondary.b3lyp.final.bgf"
g.Terminal="/qcfs/noische/research/PEI/201510/primary.b3lyp.bgf"
g.ff="/home/noische/ff/DREIDING2.21.ff"

g.num_target_node = 10
g.add_focal_point()
g.build_random()

# Add label to edges
lbl = {}
for i in g.edges():
    lbl[i] = g[i[0]][i[1]]['branch']

# Draw
#pos = nx.spring_layout(g)
#pos = nx.spectral_layout(g)
pos = nx.graphviz_layout(g, prog='twopi', args='')
nx.draw(g, pos)
nx.draw_networkx_edge_labels(g, pos, edge_labels=lbl)
plt.axis('off')
plt.show()

#g.add_monomer(0, 1)
#g.add_monomer(1, 1)
#g.add_monomer(2, 1)

#print(g.nodes())
#print(g.edges())
#print(g.calculate_wiener_index())
#g.compile("test_random_10.bgf")
