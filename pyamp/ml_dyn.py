### built-in modules
import torch
import torch.nn as nn

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

