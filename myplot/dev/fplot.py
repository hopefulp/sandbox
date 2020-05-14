#!/home/joonho/anaconda3/bin/python

import argparse
import re
from myplot2D import mplot_nvector

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

def draw_1file(onefile,ncol,ncx,title,xt,yt,Lsave):
    if title == None:
        title = re.split('\.',onefile)[0]

    with open(onefile,"r") as f:
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
    print("title {} xtitle {} ytitles {}".format(title,xt,ylabels))
    # y is plotted by column
    #x=[]
    mplot_nvector(x, y, Title=title, Xtitle=xt, Ytitle=yt, Ylabels=ylabels,Lsave=Lsave)

    return 0

def main():
    parser = argparse.ArgumentParser(description='2D plot w. list of files "-f f1 f2 f3 ..."')
    #parser.add_argument('job', choices=["ene"], help='type of input for gromacs')
    parser.add_argument('files', nargs='+',help='add all the files')
    parser.add_argument('-nc', '--ncolumn', default=1, type=int, help='number of columns in each file')
    parser.add_argument('-x', '--ncolumn_x', action='store_true', help='whether use 1st column as X')
    parser.add_argument('-t', '--title', help='title of figure would be filename')
    parser.add_argument('-xt', '--xtitle', default='t(0.2 fs)', help='X title')
    parser.add_argument('-yt', '--ytitle', default='E(eV)', help='Y title')
    parser.add_argument('-s', '--save', action='store_true', help='Save figure')
    args = parser.parse_args()
    
    if len(args.files) == 1:
        print("draw 1file with -x option and multiple y's can be drawn")
        draw_1file(args.files[0],args.ncolumn,args.ncolumn_x,args.title,args.xtitle,args.ytitle,args.save)
    else:
        draw_nfile(args.files,args.ncolumn,args.ncolumn_x,args.title,args.xtitle,args.ytitle,args.save)
    return 0

if __name__=='__main__':
	main()	


