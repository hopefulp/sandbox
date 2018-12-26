from mplt_mo_ini import *

def lprint(obj):
    if type(obj) == dict:
        for k, v in obj.items():
            print '%s : %s' % ( k, v)
            """
            if hasattr(v, '__iter__'):
                print k
                dumpclean(v)
            else:
                print '%s : %s' % ( k, v)
            """
    elif type(obj) == list:
        for v in obj:
            print v
            """
            if hasattr(v, '__iter__'):
                dumpclean(v)
            else:
                print v
            """
    else:
        print "type works for only dict and list for line print"
        print obj
    return 0        

def lprint_sorted(obj):
    if type(obj) == dict:
        a = sorted(obj.items(), key=lambda i: i[1], reverse=True)
        for k, v in a:
            if hasattr(v, '__iter__'):
                print k
                dumpclean(v)
            else:
                print '\t%-10s : %10.2f' % ( k, float(v))
    elif type(obj) == list:
        for v in obj:
            if hasattr(v, '__iter__'):
                dumpclean(v)
            else:
                print i, v
    else:
        print "type works for only dict and list for line print"
        print obj
    return 0 

def fprint_pair_list(fname, mylist):
    with open(fname, 'w') as f:
        for pair in mylist:
            f.write(str(pair[0]) + '\t' + str(pair[1]) + '\n')

def fread_pair_list(fname):
    with open(fname, 'r') as f:
        pair_list = [ [ int(x) for x in line.split()] for line in f ]
    return pair_list            
        








