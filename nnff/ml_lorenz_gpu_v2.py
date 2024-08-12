#!/home/joonho/anaconda3/bin/python

import sys
from common import whereami, list2str
import argparse
import numpy as np
import time
import os
import ga_ini
from scipy import integrate

from sklearn.model_selection import train_test_split
from torch import cuda

import torch
import torch.nn as nn
from torch import optim
import torch.nn.functional as F

Llog=0
### build MATLAB functions
class Radbas(nn.Module):
    def __init__(self, mean=0, std=1, min=0.1, max=0.9):
        super(Radbas, self).__init__()
        self.mean = mean
        self.std =std
        self.min = min
        self.max = max
        
    def forward(self, x):
        gauss = torch.exp((-(x - self.mean) ** 2)/(2* self.std ** 2))
        return torch.clamp(gauss, min=self.min, max=self.max)
        

class purelin(nn.Module):
    def __init__(self, x=0):
        super(purelin, self).__init__()
    def forward(self, x):
        return x
### Dynamics parameters
def fn_govern(fn_name):
    global dt, T, t
    dt = 0.01
    T = 8
    t = np.arange(0,T+dt,dt)
    if Llog: print(t.shape)

    if fn_name == 'lorenz':
        global beta, sigma, rho
        beta = 8/3
        sigma = 10
        rho = 28
    return 0

def lorenz_eq(x_y_z, t0, sigma, beta, rho):
    x, y, z = x_y_z
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]


### Build up DB: using dynamics of governing function, plot option
def make_db(fn_name, ndata, Lgraph=False):
    ### to plot Database
    if Lgraph:
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from mpl_toolkits.mplot3d import Axes3D
        rcParams.update({'font.size': 18})
        plt.rcParams['figure.figsize'] = [12, 12]
        fig,ax = plt.subplots(1,1,subplot_kw={'projection': '3d'})
    global t, dt, T
    ntraj = int(ndata/(T/dt)) 

    nn_input = np.zeros((ntraj*(len(t)-1),3))
    nn_output = np.zeros_like(nn_input)

    np.random.seed(123)
    x0 = -15 + 30 * np.random.random((ntraj, 3))
    ### execution test
    global x0_exe
    x0_exe = x0[-3:,:]
	### make Database: Lorenz
    if fn_name == 'lorenz':
        global beta, sigma, rho
        x_t = np.asarray([integrate.odeint(lorenz_eq, x0_j, t, args=(sigma,beta,rho)) for x0_j in x0])
        for j in range(ntraj):
            nn_input[j*(len(t)-1):(j+1)*(len(t)-1),:] = x_t[j,:-1,:]
            nn_output[j*(len(t)-1):(j+1)*(len(t)-1),:] = x_t[j,1:,:]
            if Lgraph:
                x, y, z = x_t[j,:,:].T
                ax.plot(x, y, z,linewidth=1)
                ax.scatter(x0[j,0],x0[j,1],x0[j,2],color='r')
        if Lgraph:
            ax.view_init(18, -113)
            plt.show()
    else:
        print("add more function in source")
    print(f"ntraj {ntraj}")
    return nn_input, nn_output

### Data from the above
def divide_db(nn_input, nn_output, db_ratio):
    if Llog: print(nn_input.shape, nn_output.shape)
	### Divide training, validation and test set
    ndata     = len(nn_input)
    tr_index    =     0
    ndata_tr    = int(ndata * db_ratio[0])
    ndata_val   = int(ndata * db_ratio[1])
    
    ### training data
    nn_tr_in   = nn_input[:ndata_tr]
    nn_tr_out  = nn_output[:ndata_tr]

    ### validation data
    if len(db_ratio) == 3:
        val_end = ndata_tr + ndata_val
        nn_val_in  = nn_input[ndata_tr:val_end]
        nn_val_out = nn_output[ndata_tr:val_end]
    else:
        val_end = ndata_tr
        nn_val_in  = []
        nn_val_out = []
    ### test data
    nn_te_in   = nn_input[val_end:]
    nn_te_out  = nn_output[val_end:]

    print(f"total data: {ndata}, training data: {len(nn_tr_in)}, test data: {len(nn_te_in)}")
    if len(db_ratio) == 3:
        print(f"val data: {len(nn_val_in)}")
    nn_tr=(nn_tr_in, nn_tr_out)
    nn_te=(nn_te_in, nn_te_out)
    nn_val=(nn_val_in, nn_val_out)
    return nn_tr, nn_te, nn_val

### Build Neural Network Model
def make_model(nhl):
    model = nn.Linear(3,nhl[0])
    model = nn.Sigmoid(model)
    for i in len(nhl):
        model = nn.Linear(nhl[i], nhl[i+1])(model)
        model = nn.Sigmoid(model)
    model = nn.Linear(nhl[-1],3)(model)
    model = ReLU(model)

    '''
    model = nn.Sequential(
		nn.Linear(3,nhl[0]),
		nn.Sigmoid(),
		nn.Linear(nhl[0],nhl[1]),
		Radbas(),
		nn.Linear(nhl[1],nhl[2]),
		purelin(),
		nn.Linear(nhl[2],3)
		)
    '''
    return model

class MyModel(nn.Module):
    def __init__(self, finp, fout, nnodes, act=nn.Tanh):
        super(MyModel, self).__init__()
        self.act =  act()
        self.linears = nn.ModuleList([nn.Linear(finp, nnodes[0])])
        self.linears.extend([nn.Linear(nnodes[i], nnodes[i+1]) for i in range(0, len(nnodes)-1)])
        self.out = nn.Linear(nnodes[-1], fout)

    def forward(self, x):
        for l in self.linears:
            x = self.act(l(x))
        x = self.out(x)
        return x


### Using NN Models: train, test, 
def train_model(model, X, Y, max_iter):
    optimizer = optim.Adam(model.parameters())

    for iepoch in range(max_iter):
        optimizer.zero_grad()
        hypothesis = model(X)
        loss = nn.MSELoss()
        cost = loss(hypothesis,Y)
        cost.backward()             # backward propagation
        optimizer.step()            # 
        if iepoch % 100 == 0:
            print(iepoch, cost.item())
    return model

def test_model(model, nn_data, data_kind):
    ''' 
    test model (model(X)) with X=nn_data[0] and Y=nn_data[1]
        error is average for each data (xi,yi,zi)-->(xf,yf,zf)
    '''
    loss = nn.MSELoss()
    with torch.no_grad():
        model.eval()
        X = torch.FloatTensor(nn_data[0])
        Y = torch.FloatTensor(nn_data[1])
        if cuda.is_available():
            X = X.cuda()
            Y = Y.cuda()
        y_pred = model(X)
        cost = loss(y_pred,Y)
    print(f" {data_kind} loss = {cost}")
    
    return torch.Tensor.cpu(cost)

### Running dynamics using NN Model for trajectory: also include Governing function dynamics
def run_NNdyn(model, fn_name, Lgraph):
    global t                    # defined at fn=lorenz in make_db 

    np.random.seed(100)   # original 139
    num_traj = 3

    nn_flow = np.zeros((num_traj, len(t), 3))   # flow to column
    #nn_flow[:, 0, :] = -15 + 30 * np.random.random((num_traj, 3))  # supply num_traj(3) by 3 random values
    global x0_exe
    nn_flow[:, 0, :] = x0_exe # test data comparison in loss function
    if fn_name == 'lorenz':
        global beta, sigma, rho     # defined at make_db
        x_t = np.array([integrate.odeint(lorenz_eq, nn_flow[i, 0, :], t, args=(sigma,beta, rho)) for i in range(num_traj)])

    with torch.no_grad():
        for jj, tval in enumerate(t[:-1]):
            inp = torch.Tensor(nn_flow[:,jj,:])     # becomes input of the next step-->error accumulates
            if cuda.is_available():
                inp = inp.cuda()
                nn_flow[:, jj+1, :] = torch.Tensor.cpu(model(inp))
            else:
                nn_flow[:, jj+1, :] = model(inp)     # the output of nn_flow[:, jj+1, :]
    ### 
    X = torch.FloatTensor(x_t[:,:-1,:])
    Y = torch.FloatTensor(x_t[:,1:,:])
    test_model(model, (X, Y), data_kind='Running NN dynamics')  # error for each data point

    if Lgraph:
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from mpl_toolkits.mplot3d import Axes3D
        rcParams.update({'font.size': 18})
        plt.rcParams['figure.figsize'] = [12, 12]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for j in range(num_traj):
            x, y, z = x_t[j, :, :].T
            xd, yd, zd = nn_flow[j, :, :].T                                  # xd for deep-approx
            ax.plot(x, y, z, linewidth=1)
            ax.plot(xd, yd, zd, '--', lw=1)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.scatter(x[0], y[0], z[0], color='r')
        
        ax.view_init(10, -13)
        plt.show()

        fig, axes = plt.subplots(3,3, figsize=(15,12))

        for j in range(num_traj):
            x, y, z = x_t[j, :, :].T
            xd, yd, zd = nn_flow[j, :, :].T
            axes[0][j].plot(t, x, t, xd)
            axes[1][j].plot(t, y, t, yd)
            axes[2][j].plot(t, z, label='direct solver')
            axes[2][j].plot(t, zd, label='NN')
        axes[0][0].set_title('Model1')
        axes[0][1].set_title('Model2')
        axes[0][2].set_title('Model3')
        axes[0][0].set_ylabel('X')
        axes[1][0].set_ylabel('Y')
        axes[2][0].set_ylabel('Z')
        axes[2][0].set_xlabel('t')
        axes[2][1].set_xlabel('t')
        axes[2][2].set_xlabel('t')
        plt.legend()
        print(x.shape) 

        plt.show()
    ### data test
    return 0


def NNdynamics(job, model, msave, fn_name, ndata, dt, nfeature, nhl, db_ratio, conv, max_iter, Lga, igen, Lgraph):
    if Lgraph:
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from mpl_toolkits.mplot3d import Axes3D
        rcParams.update({'font.size': 18})
        plt.rcParams['figure.figsize'] = [12, 12]
        
    ### return governing function control parameter: (global) t, dt, T, function arguments(lorenz: beta, rho, sigma)
    fn_govern(fn_name)
   
    ### data was not saved yet: calculate Lorenz by NA
    ### make database with Numerical Analysis
    nn_input, nn_output = make_db(fn_name, ndata, Lgraph)
    ### divide data into for train, [validate and ]test
    nn_tr, nn_te, nn_val = divide_db(nn_input, nn_output, db_ratio)     # returns tuple of ndarray
    ### digress by train or load
    if job == 'tr':
        #model = make_model(nhl)
        model = MyModel(3, 3, nhl)
        X = torch.FloatTensor(nn_tr[0])
        Y = torch.FloatTensor(nn_tr[1])
        ### GPU
        if cuda.is_available():
            #device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
            #print(device)
            model = model.cuda()
            X = X.cuda()
            Y = Y.cuda()
        ### train model with data
        start = time.time()
        model = train_model(model, X, Y, max_iter)
        print("time :", time.time()-start)
        ### save
        hl2str = list2str(nhl)
        if Lga:
            pre='ch'    # chromosome
            msave = f'{pre}{igen:02d}_'+hl2str+'.pt'
            ### not saving model at the moment
        else:
            if not msave:
            ### make filename for saving model
                pre='tmp'
                i=0
                msave = f'{pre}{i:02d}_'+hl2str+'.pt'
                while os.path.isfile(msave):
                    i+=1
                    msave = f'{pre}{i:02d}_'+hl2str+'.pt'
            print(f"save model {msave}")
            torch.save(model, msave)
            
    #### if job == test, Load model
    elif job == 'te' or job == 'dyn':
        if not model:
            print("load model for test: -m a.pt")
            sys.exit(10)
        else:
            model = torch.load(model)
    ### 
    ### validation and test
    if job == 'tr' or job == 'te':
        if len(db_ratio) == 3:      # do not use nn_val=([],[]) is not False
            test_model(model, nn_val, 'validate')
        if Llog: print(nn_te[0].shape)
        cost = test_model(model, nn_te, 'test')
        if Lga:
            score = ga_ini.score(cost, ref=1.0e-6)
            with open('score.dat', 'w') as f:
                f.write(' '.join("%2s" % str(nodes) for nodes in nhl))
                f.write(f" {score:10.5f}\n")
            os.system(f"cat {'score.dat'} >> {ga_ini.fitness_onegen}")

    ### execute NN-dynamics
    run_NNdyn(model, fn_name, Lgraph)
    return 0

def main():
    parser = argparse.ArgumentParser(description='Machine Learning for Lorenz equation')
    parser.add_argument('job', choices=['tr','te','dyn'],   help='name of dynamical function')
    parser.add_argument('-m', '--model', help='to load pytorch model, input model file')
    parser.add_argument('-ms', '--savefile', help='filename to save pt file')
    parser.add_argument('-fn', '--function_name',   default='lorenz',   help='name of dynamical function')
    parser.add_argument('-nd', '--ndata', type=int, default=80000,      help='total number of data')
    parser.add_argument('-dt', '--time_step', type=float, default=0.01, help='time step')
    parser.add_argument('-nx', '--num_feature', default=3, type=int, help='number of feature vector')
    parser.add_argument('-hl', '--hidden_layers', nargs='*', type=int, default=[10,10,10], help='Hidden Layer of lists of integer')
    parser.add_argument('-dbp', '--db_partition', default=3, type=int, choices=[2,3], help='part data to training:[validation:]test')
    parser.add_argument('-acc', '--accuracy', default=0.001, type=float, help='training accuracy')
    parser.add_argument('-mi', '--max_iteration', default=100000, type=int, help='stop at max_iteration for comparison')
    parser.add_argument('-ga', '--gen_alg', action='store_true', help='run genetic algorithm')  
    parser.add_argument('-ig', '--igen', default=0, type=int, help='input generation number')  
    parser.add_argument('-g', '--Lgraph', choices=['1','2','12','d','x','dx'], help='plot 1,d: database, 2,x: execution')
    args = parser.parse_args()
    
    if args.db_partition == 2:
        db_ratio=[0.85, 0.15]           # training (+validation), test
    elif args.db_partition == 3:    
        db_ratio=[0.7, 0.15, 0.15]      # training, validation, test

    NNdynamics(args.job, args.model, args.savefile, args.function_name, args.ndata, args.time_step, args.num_feature, args.hidden_layers, db_ratio, args.accuracy, args.max_iteration, args.gen_alg, args.igen, args.Lgraph)
    return 0

if __name__ == '__main__':
    main()

