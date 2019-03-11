#!/home/joonho/anaconda3/bin/python

import os

dname_full = os.path.dirname(__file__)
dname = dname_full.split('/')[-1]

comment = "This displays module and executable script in %s" % dname

print(comment)

modules = []
scripts = []

for f in os.listdir(dname_full):
    
    ext = f.split('.')[-1]
    if ext == 'py':
        #print(f)
        f_full = dname_full + '/' + f
        if os.access(f_full, os.X_OK):
            scripts.append(f)
        else:
            modules.append(f)
print("modules: ", modules)
print("scripts: ", scripts)

#for f in modules:
#    import
#    print(f.__doc__)
#for f in scripts:
#    print(f.__doc__)

