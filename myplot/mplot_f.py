#!/home/joonho/anaconda3/bin/python
''' read file and plot '''

import argparse
import re
import sys
import numpy as np
from mplot2D import mplot_nvector, mplot_twinx
from plot_job import get_jobtitle
from plot_job import Ni_4x as Nix
from my_chem import *
from common import *
import parsing 

def draw_nfile(files,ncol,title,xt,yt,Lsave):
    nfile = len(files)
    xl  = []      # x for time
    y1 = []     # y for energy
    y2 = []

    #print("%d %d" % (nfile, ncol))
    # as for f1(x,y), f2(x,y), plot x, f1(y), f2(y)
    if nfile == 2 and ncol == 2:
        with open(files[0],"r") as f, open(files[1],"r") as g:
            for line1, line2 in zip(f, g):
                line1.strip()
                line2.strip()
                l_line1 = line1.split()
                l_line2 = line2.split()
                xl.append(float(l_line1[0]))
                y1.append(float(l_line1[1]))
                y2.append(float(l_line2[1]))
    #for x, y, z in zip(xl, y1, y2):
    #    print("%10.5f%15.5f%15.5f" % (x, y, z))

    _mplot_2f3c(xl, y1, y2, files[0], files[1],title,title,xt,yt, Lsave)

    return 0

def draw_1file(file_,ncx,ncy,title,xt,yt,Lsave):
    if title == None:
        title = re.split('\.',file_)[0]

    with open(file_,"r") as f:
        lines=f.readlines()
        x = []     
        y = []     # y as 2d [ [y1], [y2], ...] 
        i = 0
        for line in lines:
            items = line.strip().split()
            if i == 0 and re.search("\w",line): # treat 1st line for title whether there is word
                if xt == None:
                    xt = items[0]
                ylabels =items[1:]              # multiple y titles
            else:
                x.append(int(items[0]))
                y_line = []     # convert string to float for a list
                for y_ in items[1:]:
                    y_line.append(float(y_))    # make 1D array
                y.append(y_line)                # make 2D array
            i+=1
    #print("title {} xtitle {} ytitles {}".format(title,xt,ylabels))
    # y is plotted by column
    #x=[]
    mplot_nvector(x, y, title, xt, yt, ylabels,Lsave)

    return 0

def convert2value(st):
    try:
        val = float(st)
    except:
        if re.search('j', st):
            val = j2cal
            if re.search('-', st):
                val *= -1
        else:
            print("Error:: No transform unit for y-scale")
            sys.exit(2)
    #print(f"return scale {val}")
    return val


def get_yscale_value(yscale):
    ''' treat yscale, values as list '''
    values=[]
    for ys in yscale:
        val = convert2value(ys)
        values.append(val)
    return values
def draw_twinx(fs,icx,x_val,icy,job,title,xt,yt,yl,Lsave,yscale,colors):
    """
    not used now
    This works for several files with only one y-value
    if x is included in file, modify
    if 1st line is label, modify
    """
    if title:
        tag_title=1
    else:
        job_title=get_jobtitle(job,title, xt, yt)
        tag_title=0

    if len(fs) != 1:
        if yl:
            if len(fs) != len(yl):
                print("Input y-label accurately")
                sys.exit(1)
        else:
            yl=[]
            for f in fs:
                pre, _ = fname_decom(f)
                yl.append(pre)

    y_2d=[]
    for f1 in fs:
        with open(f1,"r") as f:
            lines=f.readlines()
            x = []     
            y_val = []     # y as 2d [ [y1], [y2], ...] 
            i = 0
            for line in lines:
                items = line.strip().split()
                if items:
                
                    if icx:
                        x.append(float(items[icx-1]))
                    #y_line = []     # convert string to float for a list
                    if icy:
                        pass
                    ### one y-value for a file
                    else:
                    #for y in icy:
                        #y_line.append(float(items[y-1]))    # make 1D array
                        #print(f"{y} {items[y-1]}")
                        #y_val.append(float(items[y-1]))    # make 1D array
                        y_val.append(float(items[0]))    # make 1D array
                        #y.append(y_line)                # make 2D array
                i+=1
            y_2d.append(y_val)
    if not x:
        if x_val:
            x = x_val
        else:
            x=Nix
    #x=Nix
    ### y.ndim can be 1 or 2
    print(f'length of x: {len(x)}')
    if len(fs) == 1:
        y = y_val
    else:
        y = y_2d
    #yscale = get_yv_scale(yscale)           # get_yv_scale returns list
    yscale=[0.01,0.01,1.0]
    y=np.array(y)*np.array(yscale)[:,None]
    print("before draw mplot_twinx")
    print(f"{x} ")
    if tag_title:
        print(f"title = {title} xt = {xt}")
        mplot_twinx(x, y, Title=title, Xtitle=xt, Ytitle=yt, Lsave=Lsave,Ylabels=yl, Colors=colors)
    else:
        mplot_twinx(x, y, Title=job_title.title, Xtitle=job_title.xtitle, Ytitle=job_title.ytitle, Lsave=Lsave,Ylabels=yl,Colors=colors)
    return 0
def get_title(name):
    title = f_root(name)
    return title.upper()

"""
    This works for several files with only one y-value and
    one file with several y-values
    if x is included in file, modify
    if 1st line is label, modify
"""    
def draw_file(flist,icx,x_ext,icy,job,plot_dict,Lsave,yscale,twin_dict):
    '''
    flist   file list
    icy     list
    
    if plot_dict.get("xlabel"): title = plot_dict['xlabel']
    if args.ylabel: plot_dict['ylabel'] = args.ylabel
    if args.xlim:   plot_dict['xlim']   = args.xlim
    if args.title:  plot_dict['title']  = args.title
    if args.colors: plot_dict['colors'] = args.colors
    if args.legends:plot_dict['legend'] = args.legends
    twin_dict['Ltx'] = args.twinx
    if args.icfy2:  twin_dict['ic_y2']  = args.icfy2
    if args.ylabel2:twin_dict['ylabel2']= args.ylabel2  

    '''
    if not x_ext:
        x=[]
    else:
        x = x_ext
    nysize = len(x)
    ### modify get title
    if not plot_dict.get("title"):
        plot_dict['title'] = get_title(flist[0])
        
    #    tag_title=1
    #else:
    #    job_title=get_jobtitle(job,title, xlabel, ylabel)
    #    tag_title=0
    ### if 1 file, prepare legend=[]
    ### try nfile

    ### obtain legend
    if not plot_dict.get("legend"):
        legend=[]                   # each column may have legend and obtained later
        if len(flist) != 1:
            if len(flist) != len(legend):
                print("Input y-label accurately")
                sys.exit(1)
            for f in flist:
                dirname, pre, _ = dirfname_decom(f)
                legend.append(dirname)
            plot_dict['legend'] = legend
    ### in function: draw_f
    x_2d=[]
    y_2d=[]
    ### scan file list
    for f1 in flist:
        with open(f1,"r") as f:
            lines=f.readlines()
            x_val = []     
            y_val = []     # y as 2d [ [y1], [y2], ...] 
            i = 0
            ### if 1 file, get line-labels here
            for line in lines:
                items = line.strip().split()
                if len(flist) == 1 and i == 0:
                    ### if there is character, read labels
                    if parsing.is_there_char(line.strip()):
                        ncol = len(items)
                        if not icx:
                            icx=0
                        print(f"use icx {icx}")
                        x_title = items.pop(icx)
                        if not plot_dict.get("xlabel"):
                            xlabel = x_title
                        if not legend:
                            for item in items:
                                legend.append(item)
                        if not icy:
                            icy=[n for n in range(1,ncol)]
                            print(f"use icy {icy}")
                        ### in case twinx
                        if twin_dict['Ltx']:
                            ### reduce icy2 1
                            if twin_dict.get('ic_y2'):
                                icy2 = [ y2-1 for y2 in twin_dict.get('ic_y2') if icx < y2 ]
                            print(f"{twin_dict.get('ic_y2')}")
                        i += 1
                        continue
                    ### if 1st line is data, not legend
                    else:
                        x_val.append(items[icx])
                        y_val.append(items[icy[0]])                # for nfile, use 2nd column
                ### read data line from 2nd line  
                ### for several files --> no column title                      
                else:
                    if items:
                        if 'icx' in locals():
                            x_val.append(float(items[icx])) # what is difference between x and x_val
                        ### if one file, get ys
                        if len(flist) == 1:
                            for y in icy:
                                #y_line.append(float(items[y-1]))    # make 1D array
                                #print(f"{y} {items[y-1]}")
                                y_val.append(float(items[y]))    # make 1D array
                                #y.append(y_line)                # make 2D array
                        else:
                            y_val.append(float(items[1]))       # for nfile, use 2nd column
                i+=1
                
            ### y_2d: 1 y_val, nfiles
            ###       n y_val, 1 files
            y_2d.append(y_val)
            x_2d.append(x_val)                              # x cols could be different    
            #y_nfile.append(y_val)                    
    ### y.ndim can be 1 or 2
    ### x_2d 
    '''
    if not x:
        if x_val:
            x = x_val
        #else:
        #    x=Nix
    '''
    ### function: draw_f()
    ### try nfile
    if len(flist) == 1:
        ### if 1 y_val
        #y = y_val
        ### if n y_val
        print(f"size of y: {len(legend)}")
        ndata = len(y_2d[0])
        nrow  = ndata/ncol      # ?
        new_array=np.array(y_2d[0]).reshape(-1,len(legend))
        print(f"shape of y-data {new_array.shape}")
        y = [*zip(*new_array)]
        print(f"shape of y-data {np.array(y).shape}")
    ### in case several files: icy2(y-right axis) == file index of icy_right
    else:
        x = x_2d
        y = y_2d
        if twin_dict['Ltx']:
            icy2 = twin_dict.get('ic_y2')
    ### scaling with 0-th axis: y shape was ready for [y1, y2, ...,yn]
    ### scaling canbe *|+ for view in plot
    ### BE & SCF_TOTAL ['-1', '-j2cal'], NBO & BE ['-1'], 4f for EDA ['-j2cal'], NBO & CT ['-1','-j2c']
    #yscale = ['-1', '-j2cal' ]
    #yscale = [ '-j2cal' ]
    yscale = get_yscale_value(yscale)
    #yscale = [4158.65, 4158.65, 0]
    #y=np.array(y)+np.array(yscale)[:,None]
    #yscale=[0.01,0.01,1.0]
    print(f"x: {np.array(x).shape}")
    print(f"y: {np.array(y).shape} {np.array(yscale).shape} yscale={yscale}")
    if len(yscale) == 1:
        y=np.array(y)*yscale[0]
    else:
        y=np.array(y)*np.array(yscale)[:,None]
        
    if twin_dict['Ltx']:
        print("before draw mplot_twinx")
        print(f"title = {plot_dict['title']} xlabel = {plot_dict['xlabel']}")
        #yl2=[]

        #yl2.append(ylabel)
        #yl2.append(ylabel_r)
        #print(f"x {np.array(x).shape}, y {np.array(y).shape}, y2 {np.array(icy2).shape}")
        mplot_twinx(x, y, plot_dict, twin_dict, Lsave=Lsave)
    else:
        print("before draw mplot_nvector")
        print(f"title = {plot_dict.get('title')} xlabel = {plot_dict.get('xlabel')}")
        mplot_nvector(x, y, plot_dict, Lsave=Lsave)
    return 0

def draw_v(y_val, job, title, xtitle, ytitle, Lsave,yscale):
    job_title=get_jobtitle(job,title, xtitle, ytitle)

    #if job == 'eda':
    #    if yscale:
    #        ys = np.array(yscale) * j2cal
    #x=Nix
    #ys = 1
    #y=np.array(y_val)*ys
    #mplot_nvector(x, y, Title=job_title.title, Xtitle=job_title.xtitle, Ytitle=job_title.ytitle, Lsave=Lsave)
    mplot_nvector(y_val, y_val, title=job_title.title, xlabel=job_title.xtitle, ylabel=job_title.ytitle, Lsave=Lsave)
    
    return 0

def main():
    parser = argparse.ArgumentParser(description='Drawing files with values/files')

    inp = parser.add_mutually_exclusive_group()
    inp.add_argument('-v', '--values', nargs='*', help=' input y-values as y')
    inp.add_argument('-f', '--files', nargs='+',help='add all the files')
    #parser.add_argument('-xv', '--xvalues', nargs='*', help='X values')
    #g_value=parser.add_argument_group('Values', description="get argument as y-values")
    #g_value.add_argument('-yv', nargs='+', type=float, help='list y-values')
    xvalues=parser.add_mutually_exclusive_group()
    xvalues.add_argument('-icx', '--icolumn_x', type=int, default=0, help='column index of X')
    xvalues.add_argument('-x', '--x_ext', help='x-coord is supplied externally')
    g_file=parser.add_argument_group('Files', description="get input files")
    #g_file.add_argument('-nc', '--ncolumn', default=1, type=int, help='number of columns in each file')
    g_file.add_argument('-icy', '--icolumn_y', default=[1], nargs="*", type=int, help='column index of Y')
    g_file.add_argument('-ys', '--y_scale', default=['1.0'], nargs="+", help='scale factor for Y -kj2kcal')
    #g_file.add_argument('-ys', '--y_scale', nargs="*", help='scale factor for Y [value|str|str-], use for str- for "-"')
    g_twin = parser.add_argument_group('Twin-X', description='to plot using two y-axes')
    g_twin.add_argument('-tx', '--twinx', action="store_true", help='using two y-axes with twin x ticks')
    #g_twin.add_argument('-icy2', '--second_iy', default=[2], nargs="+", type=int, help='designate the index of y for 2nd y-axis')
    g_twin.add_argument('-icfy2', '--second_iy', default=[1], nargs='+', type=int, help='designate the index of y[columns or files] for 2nd y-axis')
    g_twin.add_argument('-yl2', '--ylabel2', default='E(ev)', help='input left y-axis title')
    parser.add_argument('-j', '--job', help='job of qcmo|ai|gromacs')
    plot = parser.add_argument_group(title='PLOT')
    plot.add_argument('-t', '--title', help='title of figure would be filename')
    plot.add_argument('-xl', '--xlabel', help='X title, label in mpl')
    plot.add_argument('-yl', '--ylabel', help='Y title, label in mpl')
    plot.add_argument('-xi', '--xlim', nargs=2, type=float, help='xrange xmin, xmax')
    plot.add_argument('-yi', '--ylim', nargs=2, type=float, help='yrange ymin, ymax')
    plot.add_argument('-lg', '--legend', nargs='*', help='Y labels for legend')
    plot.add_argument('-c', '--colors', nargs='*', help='Y label for legend')
    plot.add_argument('-s', '--save', action='store_true', help='Save figure')
    parser.add_argument('-u', '--usage', action='store_true', help='prints usage and exit')
    args = parser.parse_args()

    if args.usage:
        print(f"Plot from several files with different x-values\
              \n\tmplot_f.py -f Tea2SnI4spdos/TDOSF0.dat Tea2SnI4VInmspdos/TDOSF0.dat TEA2SnI4pEtOaspdos/TDOSF0.dat\
 -t TDOS -xl 'E - E$_F$ [eV]' -lg pristine V$_I$ EtO -xi -3.0 3.0 -yi 0.0 150 -c r y m\
              \n\t\
              ")
        sys.exit(0)

    ### use input values as y[]
    print(f"yscales {args.y_scale}")

    plot_dict={}
    if args.xlabel: plot_dict['xlabel'] = args.xlabel
    if args.ylabel: plot_dict['ylabel'] = args.ylabel
    if args.xlim:   plot_dict['xlim']   = args.xlim
    if args.ylim:   plot_dict['ylim']   = args.ylim
    if args.title:  plot_dict['title']  = args.title
    if args.colors: plot_dict['colors'] = args.colors
    if args.legend: plot_dict['legend'] = args.legend

    twin_dict={}
    twin_dict['Ltx'] = args.twinx
    if args.second_iy:  twin_dict['ic_y2']  = args.second_iy
    if args.ylabel2:    twin_dict['ylabel2']= args.ylabel2  

    if args.values:
        draw_v(args.values,args.job,plot_dict, args.save,args.y_scale)
    ### use input files
    elif args.files:
        print(f"read files {args.files}")
        ### x values
        if args.x_ext:
            print(f"args.x_ext = {args.x_ext}")
            print("use x-coordinate supplied by external input")
            x=Nix
        else:
            x=None
        ### n columns in 1 file, twinx: 1 column for nfiles, working ?
        print(args.second_iy)
        print(f"x before draw_file: {x}")
        draw_file(args.files, args.icolumn_x, x, args.icolumn_y,args.job,plot_dict,args.save,args.y_scale, twin_dict)
            #draw_1file(args.files[0],args.icolumn_x,args.icolumn_y,args.title,args.xlabel,args.ylabel,args.save)
    else:
        print("Error:: Turn on -v for values or -f files")
        
        #print("draw 1file with -x option and multiple y's can be drawn")
    return 0

if __name__=='__main__':
	main()	


