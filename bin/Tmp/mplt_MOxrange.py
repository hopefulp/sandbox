from mplt_qcdraw import XMIN, XMAX

#print "xrange: ", XMIN, XMAX
Xrange=XMAX-XMIN

def Xrange_1f(degen, beta):
    npart=2*degen+1
    dx=float(Xrange)/float(npart)   # devide Xrange into parts
    if beta:
        pass
    else:
        x1=XMIN+dx
        x2=x1+dx
        x3=x2+dx
        x4=x3+dx
        x5=x4+dx
        x6=x5+dx
        if degen == 1:
            return x1, x2
        elif degen == 2:
            return x1, x2, x3, x4
        elif degen==3:
            return x1, x2, x3, x4, x5, x6

    return 0

def Xrange_nf(nfile,ifile,beta,ab, degen):
    #### xrange for files
    xrange=float(XMAX-XMIN)/nfile
    #### devide xrange into partsf for degeneracy
    npart=2*degen+1
    dx=float(xrange)/float(npart)   # devide Xrange into parts

    #### xmin, xmax for a file
    if nfile==1:
        xmin=XMIN
        xmax=XMAX
    elif nfile==2:
        if ifile==0:
            xmin=XMIN
            xmax=xmin+float(xrange)
        elif ifile==1:
            xmin=XMIN+float(xrange)
            xmax=XMAX
    elif nfile==3:
        if ifile==0:
            xmin=XMIN
            xmax=xmin+float(xrange)
        elif ifile==1:
            xmin=XMIN+float(xrange)
            xmax=xmin+float(xrange)
        elif ifile==2:
            xmin=XMIN+float(xrange)*2
            xmax=xmin+float(xrange)
    #### xmin, xmax for alpha or beta in the file            
    if beta:
        x_2=float(xmax-xmin)/2
        if ab==0:
            #xmin=xmin
            xmax=xmin+x_2
        elif ab==1:
            xmin+=x_2
            #xmax=xmax
        else:
            print "error in function draw"
            exit(97)
        dx=float(x_2)/float(npart)

    x1=xmin+dx
    x2=x1+dx
    x3=x2+dx
    x4=x3+dx
    x5=x4+dx
    x6=x5+dx
    x7=x6+dx
    x8=x7+dx
    if degen == 1:
        return x1, x2
    elif degen == 2:
        return x1, x2, x3, x4
    elif degen==3:
        return x1, x2, x3, x4, x5, x6
    elif degen==4:
        return x1, x2, x3, x4, x5, x6, x7, x8
        
    return 0
    #### in the case of beta            
        
