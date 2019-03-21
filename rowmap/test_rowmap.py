from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

try:
    import _rowmap as rowmap
except ImportError:
    from build_rowmap import build_rowmap_module
    build_rowmap_module()

    # try import again
    try:
        import _rowmap as rowmap
    except:
        raise


def f(t, u):
    x1 = 0.
    xu = 1.

    n = u.size
    dx2 = ((xu - x1)/(n-1))**2
    return dx2



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


def exit_code(nexit):
    if idid == 1:
        return "computation sucessful"
    elif idid == 2:
        return "computation interrupted by solout"
    elif idid == -1:
        return "stepsize too small, i.e. hs < 10*uround*dabs(t)"
    elif idid == -2:
        return "more than iwork(1) steps"
    elif idid == -3:
        return "input is not consistent"
    elif idid == -4:
        return "internal Krylov matrix is repeatedly singular"
    else:
        return "Unknown exit code!"


def exit_details(iwork, t, hs):
    numberSteps             = iwork[4]
    numberRejected          = iwork[5]
    numberFuncEvals         = iwork[6]
    numberJacVecProds       = iwork[7]
    minLengthArrayWork      = iwork[8]
    minLengthArrayiWork     = iwork[9]

    print('NumberOfSteps: %d\nNumberRejected: %d\nNumberFuncEvals: %d\nNumberJacVecProds: %d\nMinLengthArrayWork: %d\nMinLengthArrayiWork: %d.' % (numberSteps, numberRejected,
                                        numberFuncEvals, numberJacVecProds,
                                        minLengthArrayWork,
                                        minLengthArrayiWork))
    print('Time: %.2g' % t)
    print('Step: %.2g' % hs)


def get_work(n, step_bounds = [0.25, 2.], step_safety = 0.8,
             eps = np.finfo(float).eps, ktol = 1.e-1,
             max_iter = 1000, integration_method = 2, mx = 70, lun=6):
    # prepare iwork
    iwork = np.zeros(mx + 20, np.int32)
    assert (integration_method>=1) and (integration_method<=6), ''
    iwork[0] = max_iter
    iwork[1] = integration_method
    iwork[2] = mx
    iwork[3] = lun

    # prepare work
    work = np.zeros(10 + n*(mx+11)+mx*(mx+4))
    work[0] = step_bounds[0]
    work[1] = step_bounds[1]
    work[2] = step_safety
    work[3] = eps
    work[4] = ktol

    return iwork, work


def analytic(t, x):
    return np.exp(-np.pi**2/4. * t) * np.sin(np.pi/2. * 0.5 * x)


if __name__ == '__main__':
    print(rowmap.__doc__)

    # IC
    n = 21

    x = np.linspace(0., 1., n)
    u0 = np.sin((np.pi/2.)*x)

    t0 = 0.
    tf = 2.5
    dt = 0.1
    tout = np.arange(t0, tf, dt)
    nout = n
    ncall = 0

    # sets whether rhs is autonomous or non-autonomous
    ifcn = 0

    iwork, work = get_work(n)

    atol = rtol = 1.e-4
    itol = 0
    iout = 0
    rpar = 1.
    ipar = 1
    idid = 0

    # subroutine rowmap(n,f,ifcn,t,u,tend,hs,rtol,atol,itol,ijacv,ifdt,iout,
    #                   work,lwork,iwork,liwork,rpar,ipar,idid)

    #   t,u,hs,idid = rowmap(f,ifcn,t,u,tend,hs,rtol,atol,itol,ijacv,ifdt,iout,
    #                        work,iwork,rpar,ipar,n=len(u),lwork=len(work),
    #                        liwork=len(iwork),f_extra_args=())

    print('work.sz=', work.size, ' iwork.sz=', iwork.size, ' u.sz=', u0.size)
    print('work.sz=', len(work), ' iwork.sz=', len(iwork), ' u.sz=', len(u0))

    plt.figure()
    plt.plot(x, u0, label='initial')
    #plt.plot(x, analytic(t0 + dt), label='ref')

    df = pd.DataFrame()
    df_anal = pd.DataFrame()

    def callback(*args, **kwargs):
        print('CALLBACK!')


    def solout(told, tnew, uold, unew, fold, fnew, ucon, intr, rpar, ipar):
        print('SOLOUT CALLBACK!')


    t = t0
    while True:
        nt = t + dt
        t, u, hs, iwork, idid = rowmap.rowmap(p, t, u0, t + dt, dt, rtol, atol,
                                              work, iwork, rpar, ipar, ifcn=1,
                                              iout=1,
                                              jacv=callback, fdt=callback,
                                              solout=solout)

        df[t] = u.copy()
        df_anal[t] = analytic(nt, x)
        error = np.linalg.norm(df[t] - df_anal[t], ord=2)

        print('Error at %.2g is %.2g.' % (nt, error))

        if idid is not 1 or t > tf:
            break

    print('Completed integration: %s.' % exit_code(idid))
    exit_details(iwork, t, hs)

    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(121, projection='3d')
    ax.set_title('Numerical')

    for index, row in df.iterrows():
        ax.plot(tout, x[index] * np.ones_like(row.values), row.values, color='k')

    ax.set_xlabel('Time')
    ax.set_ylabel('Space')
    ax.set_zlabel('u(x,t)')

    ax = fig.add_subplot(122, projection='3d')
    ax.set_title('Analytical')

    for index, row in df_anal.iterrows():
        ax.plot(tout, x[index] * np.ones_like(row.values), row.values, color='k')

    ax.set_xlabel('Time')
    ax.set_ylabel('Space')
    ax.set_zlabel('u(x,t)')


    plt.show()



