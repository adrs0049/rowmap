import numpy as np
from mol.MOL import MOL
from tdr.TDR import TDR

if __name__ == '__main__':
    # setup domain
    L = 10
    cellsPerUnitLength = 2**5
    h = 1. / (2**5)
    N = L * cellsPerUnitLength
    x = np.linspace(0, L, N)
    y1 = np.sin(2. * np.pi * x / L)
    y2 = np.cos(2. * np.pi * x / L)
    y0 = np.vstack((y1, y2))
    print('y0=',y0.shape)

    nop = 1

    ngb = [{ 'left' : 1, 'right' : 1 }]
    dX  = np.array([[h]])
    x0  = np.array([[0]])
    N   = np.array([[L * cellsPerUnitLength]]).astype(int)

    # number of equations
    size = 2

    # transition matrices
    trans    = np.eye(2)
    #np.zeros(size * size).reshape(size, size)
    Adhtrans = np.zeros(size * size).reshape(size, size)

    # times
    t0 = 0.
    tf = 1.
    dt = 0.1

    solver = MOL(TDR,
                 y0,
                 nop = nop,
                 ngb = ngb,
                 dX  = dX,
                 x0  = x0,
                 N   = N,
                 noPDEs = size,
                 transitionMatrix = trans,
                 AdhesionTransitionMatrix = Adhtrans,
                 t0 = t0,
                 tf = tf,
                 dt = dt
                )

    solver.run()


