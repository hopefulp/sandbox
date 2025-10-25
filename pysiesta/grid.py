import numpy as np

class Grid(object):

    def __init__(self, cell, mesh, grid):

        self._cell = cell
        self._mesh = mesh
        self._grid = grid

        nsm = mesh[0] * mesh[1] * mesh[2]
        vol = abs(np.dot(np.cross(cell[0],cell[1]), cell[2]))

        self._vol = vol
        self._dvol = vol / nsm

        dvector = np.copy(cell)
        dvector[0] = cell[0] / mesh[0]
        dvector[1] = cell[1] / mesh[1]
        dvector[2] = cell[2] / mesh[2]
        self._dvector = dvector

    def __add__(self, other):

        if isinstance(other, Grid):
            if ((self._cell == other._cell).all() and
                (self._mesh == other._mesh).all()):

                cell = self._cell
                mesh = self._mesh

                grid1 = self._grid
                grid2 = other._grid
                grid3 = grid1 + grid2

                return Grid(cell, mesh, grid3)

    def __sub__(self, other):

        if isinstance(other, Grid):
            if ((self._cell == other._cell).all() and
                (self._mesh == other._mesh).all()):

                cell = self._cell
                mesh = self._mesh
        
                grid1 = self._grid
                grid2 = other._grid
                grid3 = grid1 - grid2

                return Grid(cell, mesh, grid3)

    def __mul__(self, other):

        if isinstance(other, Grid):
            if ((self._cell == other._cell).all() and
                (self._mesh == other._mesh).all()):

                cell = self._cell
                mesh = self._mesh
                dvol = self._dvol

                grid1 = self._grid
                grid2 = other._grid
                integ = np.sum(grid1 * grid2 * dvol)

                return integ

    def sum(self):

        dvol = self._dvol
        grid = self._grid
        return np.sum(grid * dvol)


    def pav(self):

        from scipy import interpolate

        
