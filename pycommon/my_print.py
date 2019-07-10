from mplt_mo_ini import *

def lprint(obj):
    if type(obj) == dict:
        for k, v in obj.items():
            print ('%s : %s' % ( k, v))
            """
            if hasattr(v, '__iter__'):
                print k
                dumpclean(v)
            else:
                print '%s : %s' % ( k, v)
            """
    elif type(obj) == list:
        for v in obj:
            print (v)
            """
            if hasattr(v, '__iter__'):
                dumpclean(v)
            else:
                print v
            """
    else:
        print ("type works for only dict and list for line print")
        print (obj)
    return 0        

def lineprint_sorted(obj):
    if type(obj) == dict:
        #a = sorted(obj.items(), key=lambda i: i[1], reverse=True)
        #for k, v in a:
        for k, v in obj:
            if hasattr(v, '__iter__'):
                print (k)
                lineprint_sorted(v)
            else:
                print ('\t%-10s : %10.2f' % ( k, float(v)))
    elif type(obj) == list:
        for v in obj:
            if hasattr(v, '__iter__'):
                lineprint_sorted(v)
            else:
                print (i, v)
    else:
        print ("type works for only dict and list for line print")
        print (obj)
    return 0 

def fprint_pair_list(fname, mylist):
    with open(fname, 'w') as f:
        for pair in mylist:
            f.write(str(pair[0]) + '\t' + str(pair[1]) + '\n')

def fread_pair_list(fname):
    with open(fname, 'r') as f:
        pair_list = [ [ int(x) for x in line.split()] for line in f ]
    return pair_list            

def fprint_pair_2Dlist(fname, mylist):
    with open(fname, 'w') as f:
        for links in mylist:
            for pair in links:
                f.write(str(pair[0]) + '\t' + str(pair[1]) + '\n')
            f.write("\n")
        
def fread_pair_2Dlist(fname):
    pair_2Dlists=[]
    with open(fname, 'r') as f:
        one_link=[]
        for line in f:
            pair=line.split()
            if pair:
                one_link.append(pair)
            else:
                pair_2Dlists.append(one_link)
                one_link=[]
        pair_2Dlists.append(one_link)
    return pair_2Dlists






