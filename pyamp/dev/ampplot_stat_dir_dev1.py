#!/home/joonho/anaconda3/bin/python

import argparse
import sys
import glob
import pickle
import re
from common import f_parsing
import amp_ini

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def namestr(obj, namespace):
    '''
    name in namespace should have the unique value
    '''
    #print(f"obj = {obj} with namespace {namespace}")
    return [ name for name in namespace if namespace[name] == obj ]

def amp_stat_plot(dirs, fname, legend_style, HL, E_conv,f_conv,f_coeff, ntrain, ntest, xlabel=None, ylabel=None, title=None):
    if not dirs:
        dirs = glob.glob("nd*")
    dirs.sort(key=natural_keys)
    #print(dirs)
    global mse, bias_sq, varx_y, varx, vary, cov, r_corr, rmse, maxres
    mse     = []
    bias_sq = []
    varx_y  = []
    varx    = []
    vary    = []
    cov     = []
    r_corr  = []
    rmse    = []
    maxres  = []

    x = []
    print(f"dirs = {dirs}")
    for d in dirs:
        x.append(int(d[5:]))
        with open(d+'/'+fname, 'rb') as f:
            p = pickle.load(f)
        for key in p.keys():
            #print(key)
            exec(f"{key}.append({p[key]})")
    print(f"x[] = {x}")

    ### Plot components
    if legend_style == 'mse':
        y=[mse, bias_sq, varx_y, r_corr]
        y_ax2 = [3]
    elif legend_style == 'maxres':
        y = [rmse, maxres]
        y_ax2 = []
    legend=[]
    for name in y:
        name_list = namestr(name, globals())
        #print(f"get name_list {name_list} from function namestr()")
        legend.append(name_list[0].upper())
            
    if True:
        modulename='myplot2D'
        if modulename not in sys.modules:
            import myplot2D
        print(f"ylabel {ylabel}")
        if not ylabel: 
            ylabel = ['(meV)^2', 'Correlation']
        elif len(ylabel) == 1:
            ylabel = ylabel[0]
        if not xlabel: xlabel = 'NData(tr)'
        if y_ax2:
            myplot2D.mplot_twinx(x, y, y_ax2, title = title, xlabel=xlabel, ylabel=ylabel, legend=legend)
        else:
            myplot2D.mplot_nvector(x, y, title = title, xlabel=xlabel, ylabel=ylabel, legend=legend)
    return 0

def main():
    parser = argparse.ArgumentParser(description='plot amp test using my2dplot')
    parser.add_argument('dirs', nargs='*', help='amp dirs for statistics')
    parser.add_argument('-f', '--fname',  help='energy file for amp test')
    parser.add_argument('-le', '--legend',  default='mse', choices=['mse','maxres'], help='select legend components')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[4], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_conv', default=0.0001, type=float, help='energy convergence limit')
    parser.add_argument('-fl', '--f_conv', default=0.1, type=float, help='force convergence limit')
    parser.add_argument('-fc', '--f_coeff', default=0.04, type=float, help='force coefficient')
    parser.add_argument('-ntr', '--ndata_train', type=int, help="N total data")
    parser.add_argument('-nte', '--ndata_test', type=int, help="N test data")
    parser.add_argument('-xl', '--xlabel', help='xlabel for matplot')
    parser.add_argument('-yl', '--ylabel', nargs='*', help='ylabel for matplot')
    parser.add_argument('-tl', '--title', help='title for matplot')
    args = parser.parse_args()

    if args.title and not args.fname:
        if re.match('E', args.title):
            inf = amp_ini.ampout_te_e_chk
        elif re.match('F', args.title):
            inf = 'test_fstat_acc.pkl'
        else:
            print("input -f input file")
            sys.exit(1)
    elif args.fname:
        inf = args.fname
        

    amp_stat_plot(args.dirs,inf,args.legend,args.hidden_layer,args.e_conv,args.f_conv,args.f_coeff,args.ndata_train,args.ndata_test,args.xlabel,args.ylabel,args.title)
    return 0

if __name__ == '__main__':
    main()

