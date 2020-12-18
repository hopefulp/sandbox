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

def amp_stat_plot(dirs, dprefix,subdir, fname, keyv, legend_style,xlabel, ylabel, title, HL, E_conv,f_conv,f_coeff, ntrain, ntest):
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
    imaxres = []

    if dirs: 
        pass
    elif dprefix:
        dirs = glob.glob(f"{dprefix}*")
    dirs.sort(key=natural_keys)
    x = []
    print(f"dirs = {dirs}")
    for d in dirs:
        ### make x labels
        #ndata = ''.join(filter(lambda i: i.isdigit(), d))    # how to retrieve digits from filter object
        ### obtain integer from dirname for x-values in plot
        if dprefix:
            dir_num = re.sub(dprefix, "", d)
        else:
            dir_num = re.sub("\D", "", d)
        print(f"length of dir_num {dir_num} = {len(dir_num)}")
        if keyv == 'hl':
            dir_sub = dir_num[int(len(dir_num)/2):]
        else:
            dir_sub = dir_num
        x.append(int(dir_sub))
        ### if data locate in subdirectory
        if subdir:
            dname = f"{d}/{subdir}/{fname}"
        else:
            dname = "{d}/{fname}"
        with open(dname, 'rb') as f:
            p = pickle.load(f)
        for key in p.keys():
            #print(key)
            exec(f"{key}.append({p[key]})")
    print(f"x[] = {x}")

    ### Plot components
    if legend_style == 'mse':
        if fname == amp_ini.ampout_te_e_chk:
            y=[mse, bias_sq, varx_y, r_corr]
            y_ax2 = [3]
        ### if force
        else:
            y=[mse, r_corr]
            y_ax2 = [1]
    elif legend_style == 'maxres':
        y = [rmse, maxres]
        y_ax2 = [1]
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
    parser.add_argument('-d', '--dirs', nargs='*', help='explicit derectories')
    parser.add_argument('-p', '--dprefix', help='prefix for directories')
    parser.add_argument('-sd', '--sub_dir', help='read subdirectory below the directories with prefix')
    parser.add_argument('-f', '--fname',  help='energy file for amp test')
    parser.add_argument('-k', '--keyvalue', choices=['hl'], help='input directory properties')
    amp_gr = parser.add_argument_group(title='AMP')
    amp_gr.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[4], help='Hidden Layer of lists of integer')
    amp_gr.add_argument('-el', '--e_conv', default=0.0001, type=float, help='energy convergence limit')
    amp_gr.add_argument('-fl', '--f_conv', default=0.1, type=float, help='force convergence limit')
    amp_gr.add_argument('-fc', '--f_coeff', default=0.04, type=float, help='force coefficient')
    amp_gr.add_argument('-ntr', '--ndata_train', type=int, help="N total data")
    amp_gr.add_argument('-nte', '--ndata_test', type=int, help="N test data")
    plt_gr = parser.add_argument_group(title='PLOT')
    plt_gr.add_argument('-le', '--legend',  default='mse', choices=['mse','maxres'], help='select legend components')
    plt_gr.add_argument('-xl', '--xlabel', help='xlabel for matplot')
    plt_gr.add_argument('-yl', '--ylabel', nargs='*', help='ylabel for matplot')
    plt_gr.add_argument('-tl', '--title', help='title for matplot')
    args = parser.parse_args()

    if re.match('E', args.title) or re.search('e', args.fname, re.IGNORECASE):
        inf = amp_ini.ampout_te_e_chk
    elif re.match('F', args.title) or re.search('f', args.fname, re.IGNORECASE):
        inf = amp_ini.ampout_te_f_chk
    #elif args.fname:
    #    inf = args.fname

    amp_stat_plot(args.dirs,args.dprefix,args.sub_dir,inf,args.keyvalue,args.legend,args.xlabel,args.ylabel,args.title, args.hidden_layer,args.e_conv,args.f_conv,args.f_coeff,args.ndata_train,args.ndata_test)
    return 0

if __name__ == '__main__':
    main()

