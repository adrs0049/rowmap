import numpy as np

class Parameters(object):
    def __init__(self):
        self.D = 1.

def ProbFTrans(params)
    return params.D


"""
    Space discretization is cubic essentially, values u_i are the average
    value of a lattice cell.

    Then fluxes are computed flowing between the different cubic lattice
    volumes.

    Therefore we have to interpolate the values at the boundaries of the
    lattice

"""

def tdrComputeFacedata(t, patchId):
    pass

# compute diffusion approximation on patch patchId
def tdrDiff(t, patchId):
    pass


class Patch(object):
    def __init__(self, t, y, N, n = 1, bW = 1):
        """
            N = [N_x, N_y]
            n = Number of concentration fields
            bW = boundary width

            y    = [N_x + 2 bW, N_y + 2 bW, n]
            ydot = [N_x, N_y, n]

        """

        self.t = t

        if N[1] == 1:
            ny = 1
        else:
            ny = N[1] + 2 * bW

        self.y = np.empty((N[0] + 2 * bW, ny, n))
        self.y[:] = np.NaN

        self.y[bW + np.arange(0, N[0]), :] = np.reshape(y, (n, N[0])).T

        self.ydot = np.zeros((N[0], N[1], n))

        # TODO
        self.ComputeFaceData
        self.uDx
        self.uDy
        self.uAvx
        self.uAvy

    def clear(self):
        pass


class Domain(object):
    def __init__(self, n):

        # number of patches
        self.nop = 1

        # number of conc fields
        self.n = n

        # dimensions
        self.dim = 2

        # neighbours of each patch
        self.ngp = np.zeros((self.nop, 4))

        # neighbourhood map
        left = 0
        right = 1
        bottom = 2
        top = 3

        # lower left coordinate of each patch
        self.coord = np.zeros((self.nop, self.dim))

        # grid information for patches
        # we use a uniform rectangular grid on each patch
        self.dX = np.zeros((self.nop, self.dim))

        # number of rectangles
        self.N = np.zeros((self.nop, self.dim))

        # width of boundary
        self.boundaryWidth = 1.

        # map to patches
        self.patches = {}


    def set_additional_data(is1D=True):

        assert is1D, 'only 1d is supported atm!'

        self.ps = np.array([1])
        self.pe = np.array([self.n * self.N[0]])
        self.gridsize = self.N[0]
        self.yc = self.coord[0] + self.dX[0] * (np.arange(1, self.N[0] + 1, 1) - 0.5)
        self.cellCentreMatrix = self.yc

    def prepare_patch(self, t, y, patchId):

        assert (is1D and patchId == 0), 'this doesnt work in 1D'

        self.patches.pop(patchId, None)

        ystart = self.ps[patchId]
        yend = self.pe[patchId]

        self.patches[patchId] = Patch(t, self.N, n = self.n, bW=self.boundaryWidth)
        patch = self.patches[patchId]

        ngbPatchId = self.ngb(patchId, left)
        Cb = np.arange(0, self.boundaryWidth)
        nCb = np,arange(-self.boundaryWidth, 0)
        if ngbPatchId == patchId:
            # periodic BC
            patch.y[Cb, :] = patch.y[nCb + self.N[0] + self.boundaryWidth, :]
        else:
            assert False, 'Not implemented'

        ngbPatchId = self.ngb(patchId, right)
        if ngbPatchId == patchId:
            # periodic BC
            patch.y[self.boundaryWidth + self.N[0] + Cb, :] = patch.y[self.boundaryWidth + Cb, :]
        else:
            assert False, 'Not implemented'

        # 1D stuff is DONE HERE





