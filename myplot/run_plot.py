#!/home/joonho/anaconda3/bin/python

import argparse
import re
import sys
import numpy as np
from myplot2D import mplot_nvector, mplot_twinx
#from plot_job import get_jobtitle, Ni_6x
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


def get_yv_scale(yscale):
    ''' treat yscale, values as list '''
    values=[]
    for ys in yscale:
        val = convert2value(ys)
        values.append(val)
    return values
def draw_twinx(fs,icx,x_val,icy,job,title,xt,yt,yl,Lsave,yscale,colors):
    """
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
            x=Ni_6x
    #x=Ni_6x
    ### y.ndim can be 1 or 2
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



"""
    This works for several files with only one y-value and
    one file with several y-values
    if x is included in file, modify
    if 1st line is label, modify
"""    
def draw_f(fs,icx,x_col,icy,job,title,xt,yt,yl,Lsave,yscale,colors,Ltwinx):
    if x_col:
        x=Ni_6x
    else:
        x=[]
    if title:
        tag_title=1
    else:
        job_title=get_jobtitle(job,title, xt, yt)
        tag_title=0

    if len(fs) == 1:
        if not yl:
            yl=[]
    else:
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
            x_val = []     
            y_val = []     # y as 2d [ [y1], [y2], ...] 
            i = 0
            for line in lines:
                ### if 1 file with labels, obtain yl
                if i==0:
                    ### if there is character, read labels
                    if parsing.is_there_char(line.strip()):
                        items = line.strip().split()
                        ncol = len(items)
                        if not icx:
                            icx=1
                            print(f"use icx {icx}")
                        x_title = items.pop(icx-1)
                        if not xt:
                            xt = x_title
                        if not yl:
                            for item in items:
                                yl.append(item)
                        if not icy:
                            icy=[x for x in range(2,ncol+1)]
                            print(f"use icy {icy}")
                        i += 1
                        continue
                ### read data line from file                        
                items = line.strip().split()
                if items:
                    if icx:
                        x.append(float(items[icx-1]))
                        #x_value.append(float(items[icx-1]))
                    #y_line = []     # convert string to float for a list
                    for y in icy:
                        #y_line.append(float(items[y-1]))    # make 1D array
                        #print(f"{y} {items[y-1]}")
                        y_val.append(float(items[y-1]))    # make 1D array
                        #y.append(y_line)                # make 2D array
                i+=1
            ### y_2d: 1 y_val, nfiles
            ###       n y_val, 1 files
            y_2d.append(y_val)                    
            #y_nfile.append(y_val)                    
    ### y.ndim can be 1 or 2
    if len(fs) == 1:
        ### if 1 y_val
        #y = y_val
        ### if n y_val
        print(f"size of y: {len(yl)}")
        ndata = len(y_2d[0])
        nrow  = ndata/ncol
        new_array=np.array(y_2d[0]).reshape(-1,len(yl))
        print(f"shape of y-data {new_array.shape}")
        y = [*zip(*new_array)]
        print(f"shape of y-data {np.array(y).shape}")
    else:
        y = y_2d
    ### scaling with 0-th axis: y shape was ready for [y1, y2, ...,yn]
    ### scaling canbe *|+ for view in plot
    #yscale = get_yv_scale(yscale)
    yscale = [4158.65, 4158.65, 0]
    y=np.array(y)+np.array(yscale)[:,None]
    if tag_title:
        print(f"title = {title} xt = {xt}")
        mplot_nvector(x, y, Title=title, Xtitle=xt, Ytitle=yt, Lsave=Lsave,Ylabels=yl,Colors=colors,Ltwinx=Ltwinx)
    else:
        mplot_nvector(x, y, Title=job_title.title, Xtitle=job_title.xtitle, Ytitle=job_title.ytitle, Lsave=Lsave,Ylabels=yl,Colors=colors,Ltwinx=Ltwinx)
    return 0

def draw_v(y_val, job, title, xtitle, ytitle, Lsave,yscale):
    job_title=get_jobtitle(job,title, xtitle, ytitle)

    if job == 'eda':
        if yscale:
            ys = np.array(yscale) * j2cal
    x=Ni_6x
    y=np.array(y_val)*ys
    mplot_nvector(x, y, Title=job_title.title, Xtitle=job_title.xtitle, Ytitle=job_title.ytitle, Lsave=Lsave)
    
    return 0

def main():
    parser = argparse.ArgumentParser(description='Drawing files with values, files ')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-v', '--values', action='store_true', help='use input values as y')
    group.add_argument('-f', '--files', nargs='+',help='add all the files')

    group_v=parser.add_argument_group('Values', description="get argument as y-values")
    group_v.add_argument('-y', nargs='+', type=float, help='list y-values')
    
    group_f=parser.add_argument_group('Files', description="get input files")
    #group_f.add_argument('-nc', '--ncolumn', default=1, type=int, help='number of columns in each file')
    group_f.add_argument('-icx', '--icolumn_x', type=int, help='column index of X')
    group_f.add_argument('-x', '--x_col', help='x-coord is supplied externally')
    group_f.add_argument('-icy', '--icolumn_y', nargs="*", type=int, help='column index of Y')
    group_f.add_argument('-tx', '--twinx', action="store_true", help='using two y-axes with twin x ticks')
    group_f.add_argument('-ys', '--y_scale', default=[1], nargs="+", help='scale factor for Y [value|str|str-], use for str- for "-"')
    #group_f.add_argument('-ys', '--y_scale', nargs="*", help='scale factor for Y [value|str|str-], use for str- for "-"')
    parser.add_argument('-j', '--job', help='job of qcmo|ai|gromacs')
    parser.add_argument('-t', '--title', help='title of figure would be filename')
    parser.add_argument('-xt', '--xtitle', help='X title')
    parser.add_argument('-xv', '--xvalues', nargs='*', help='X values')
    parser.add_argument('-yt', '--ytitle', help='Y title')
    parser.add_argument('-yl', '--ylabel', nargs='*', help='Y label for legend')
    parser.add_argument('-c', '--colors', nargs='*', help='Y label for legend')
    parser.add_argument('-s', '--save', action='store_true', help='Save figure')
    args = parser.parse_args()

    ### use input values as y[]
    if args.values:
        draw_v(args.y,args.job,args.title,args.xtitle,args.ytitle,args.save,args.y_scale)
    ### use input files
    elif args.files:
        ### use twin y-axis
        #if args.twinx:
        #    draw_twinx(args.files,args.icolumn_x,args.xvalues,args.icolumn_y,args.job,args.title,args.xtitle,args.ytitle,args.ylabel,args.save,args.y_scale,args.colors)
        ### user one y-axis
        #else:
        draw_f(args.files,args.icolumn_x,args.x_col,args.icolumn_y,args.job,args.title,args.xtitle,args.ytitle,args.ylabel,args.save,args.y_scale, args.colors, args.twinx)
            #draw_1file(args.files[0],args.icolumn_x,args.icolumn_y,args.title,args.xtitle,args.ytitle,args.save)
    else:
        print("Error:: Turn on -v for values or -f files")
        
        #print("draw 1file with -x option and multiple y's can be drawn")
    return 0

if __name__=='__main__':
	main()	


