### built-in modules
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

