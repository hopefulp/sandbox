#!/gpfs/home/joonho/anaconda3/bin/python

import argparse

### nodes in Y. M. Rhee's group
def skylake_list():
    a = []
    for i in range(1,17):
        st = "node{:02d}".format(i)
        a.append(st)
    return a
    
nodes_name=['skylake', 'broadwell', 'haswell', 'other', 'sandy']
nodes_skylake = skylake_list()
nodes_broadwell = ['opt02']
nodes_haswell = ['opt05', 'slet07']
nodes_other = ['opt06', 'opt08', 'opt09']
nodes_sandy = ['opt03', 'opt04', 'opt07', 'opt10', 'slet01', 'slet02', 'slet04', 'slet06']
nodes_all = {  'skylake': nodes_skylake, 'broadwell': nodes_broadwell, 'haswell':nodes_haswell, \
                'other': nodes_other, 'sandy':nodes_sandy }

def main():
    parser = argparse.ArgumentParser(description='show nodes in SGE')
    parser.add_argument('-n', '--nodename')
    args = parser.parse_args()

    if args.nodename == None:
        print("All nodes: ", nodes_all.keys())
        for keys, values in nodes_all.items():
            print(keys)
            print(values)


if __name__ == '__main__':
    main()
    
