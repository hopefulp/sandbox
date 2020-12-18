#!/home/joonho/anaconda3/bin/python

import argparse
#import myplot2D                # qsub can't read
from common import f_parsing
import sys
from amp_util import *

def get_title(fname, HL, E_conv, f_list, ntrain, ndata, title=None):
    if f_list:
        f_conv, f_coeff = decom_ef(f_list)
    if not title:
        title = fname.split(".")[0] + "\n"
    hl = '$\\times$'.join(str(x) for x in HL)
    suptitle = "\n\nAMP Model(HL={}".format(hl) + ", "
    suptitle += "E_rms={} ".format(E_conv) + " "
    if f_list:
        suptitle += "F_err={},F_coeff={}):".format(f_conv,f_coeff) + " "
    suptitle += "test/train={}/{}".format(ndata, ntrain)
    return title, suptitle

def amp_pyplot(fdata, HL, E_conv,f_list, ntrain, ntest, title):
    title, subtitle = get_title(fdata, HL, E_conv,f_list, ntrain, ntest, title)

### file parsing
    Ys = f_parsing(fdata)
### 2 columns    
    y       = Ys[0][:]
    y_bar   = Ys[1][:]

    if not ntest:
        ntest = len(y)
    import myplot2D
    ### use Ltwinx for error axis
    #err, res_max = myplot2D.draw_2subdots(y, y_bar, title, subtitle, Ltwinx=True)
    err, res_max = myplot2D.draw_amp_twinx(y, y_bar, title, subtitle, Ltwinx=True)
    return 0

def main():
    parser = argparse.ArgumentParser(description='plot amp test using my2dplot')
    parser.add_argument('-f', '--inf', default='test_energy.dat', nargs='?', help='amp test data file ')
    #parser.add_argument('-f', '--fname',  default='OUTCAR', help='amp pot out file ')
    parser.add_argument('-hl', '--hidden_layer', nargs='*', type=int, default=[8,8,8], help='Hidden Layer of lists of integer')
    parser.add_argument('-el', '--e_conv', default=0.001, type=float, help='energy convergence limit')
    #parser.add_argument('-fl', '--f_list', nargs='*', default=[0.01,0.04], type=float, help='force convergence limit & f coeff')
    parser.add_argument('-fl', '--f_list', nargs='*', type=float, help='force convergence limit & f coeff')
    parser.add_argument('-ntr', '--ndata_train', type=int, help="N total data")
    parser.add_argument('-nte', '--ndata_test', type=int, help="N test data")
    parser.add_argument('-tl', '--title', help="title")
    
    args = parser.parse_args()

     
    amp_pyplot(args.inf,args.hidden_layer,args.e_conv,args.f_list,args.ndata_train,args.ndata_test,args.title)
    return 0

if __name__ == '__main__':
    main()

