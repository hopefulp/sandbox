#!/home/joonho/anaconda/bin/python
import sys

OUTCAR = str(sys.argv[1])
EIGENVAL = str(sys.argv[2])
ym = float(sys.argv[3])
yx = float(sys.argv[4])

def ShowTnBand(OUTCAR,EIGENVAL):
    """
    
    """

    import matplotlib.pyplot as plt
    import numpy as np

    #initiation
    #Data_of_Band : [[Fermi level], [K-point Min, K-point MAX], Energy min, Energy max]
    #, [The number of obitals, The number of spin components, The number of K-points]
    #xk_points : The list of all of k-points
    #bands : the list of eigenvalue energy
    Data_of_Band = []
    xk_points = []
    bands = []

    #Import OUTCAR file to get the infomation of fermi energy
    f = open(OUTCAR)
    list_lines = []

    #Get OUTCAR information
    for line in f.readlines():
        list_lines.append(line)

    #Get Fermi-energy information in the OUTCAR file
    a = list_lines.index(' average (electrostatic) potential at core\n')
    b = list_lines.index(' k-point     1 :       0.5000    0.5000    0.0000\n')
    list_lines = list_lines[a:b]
    Data_of_Band.append([])
    for i in range(0,len(list_lines)):
        temp=list_lines[i].split()
        if temp != []:
            if temp[0] == 'E-fermi':
                Data_of_Band[0].append(float(temp[2]))
    print (Data_of_Band)

    #==========================================================================#

    #Import EIGENVAL file to get the information of kpoint, energy etc.
    f = open(EIGENVAL)
    list_line = []

    #Get EIGENVAL information
    for line in f.readlines():
        list_line.append(line.split())
    list_line.append([])
    total_data = []
    temp = []
    for line in list_line:
        temp.append(line)
        if line == []:
            total_data.append(temp[:-1])
            temp = []
    #assign each k-point to xkpoints
    for i in range(1,len(total_data)):
        xk_points.append(float(total_data[i][0][2]))
        del total_data[i][0]

    #assign each energy according to spin compoents
    spin_comp = len(total_data[1][0])
    for i in range(1, spin_comp):
        bands.append([])
        for j in range(0, len(total_data[1])):
            bands[i-1].append([])
            for k in range(1, len(total_data)):
                bands[i-1][j].append(float(total_data[k][j][i]))
    
    #==========================================================================#
    
    #Set range of x axis(y axis==>need to update)
    xmax = max(xk_points)
    xmin = min(xk_points)
    ymin = ym
    ymax = yx
    
    #Drawing structue of fiqure to each spin components
    for i in range(0,spin_comp-1):
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        #plt.title('spin components = %s' % str(i+1))
        fig1.set_ylabel('Energy [eV]', fontsize='x-large')
        fig1.set_xlabel('')
        fig1.set_xlim(xmin, xmax)
        fig1.set_ylim(ymin, ymax)
        print ('spin components = %s' % str(i+1))
    #Poltting band structure
        for k in range(0,len(bands[i])):
            list_x = xk_points
            list_B = bands[i][k]
            array_x = np.array(list_x); array_B = np.array(list_B)
            fig1.plot(array_x, array_B, linewidth = '2', color = 'black')
            del list_B; del list_x
    plt.show()

ShowTnBand(OUTCAR,EIGENVAL)