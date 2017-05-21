from scipy.integrate import ode
import rowmap_ode_runner
import numpy as np
import pandas as pd

def p(t, u):
    xu = 1.
    x1 = 0.
    n = u.size
    ut = np.zeros_like(u)
    dx2 = ((xu-x1)/(n-1))**2

    for i in range(0, n):
        if i == 0:
            ut[i] = 0.
        elif i == n-1:
            ut[i] = 2.*(u[i-1]-u[i])/dx2
        else:
            ut[i] = (u[i+1]-2.*u[i]+u[i-1])/dx2

    return ut


if __name__ == '__main__':
    # points to solve pde on
    n = 21

    x = np.linspace(0., 1., n)
    u0 = np.sin((np.pi/2.)*x)

    t0 = 0.
    tf = 2.5
    dt = 0.1
    tout = np.arange(t0, tf, dt)

    # setup ode system
    ode15 = ode(p).set_integrator('rowmap', method='grk4t')
    ode15.set_initial_value(u0, t0)

    # result storage
    df = pd.DataFrame()

    while ode15.successful() and ode15.t < tf:
        df[ode15.t] = ode15.integrate(ode15.t + dt)




