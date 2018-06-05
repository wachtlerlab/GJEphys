import numpy as np
from scipy.optimize import basinhopping
import sys
import json
from multiprocessing import Pool

assert len(sys.argv) == 2, 'Only one argument, the path of the swcfile expected, ' + str(len(sys.argv)) + 'found'

parFile = sys.argv[1]
with open(parFile, 'r') as fle:
    pars = json.load(fle)

def accept_test(f_new, x_new, f_old, x_old):
        [A1, onset1, tau1, A2, onset2, tau2, offset] = x_new

        return all([x > 0 for x in [A1, A2, tau1, tau2]]) and bool(onset1 >= -10) and bool(onset1 <= 15) \
               and bool(onset2 >= 0) and bool(onset2 <= 40) \
               and bool(offset <= 2) and bool(offset >= -2)

        # return all([x >= 0 for x in [A1, A2, tau1, tau2]])


class ArgGenIterator:

    def __init__(self, xData, yData, p0s):


        self.xData = xData
        self.yData = yData
        self.p0s = p0s
        self.pointsDone = 0

    def __iter__(self):

        self.pointsDone = 0
        return self

    def next(self):

        if self.pointsDone < len(self.p0s):
            toReturn = (self.xData, self.yData, self.p0s[self.pointsDone])
            self.pointsDone += 1
            return toReturn
        else:
            raise StopIteration

def fitAlpha(y):

    xData, yData, p0 = y
    xData = np.array(xData)
    yData = np.array(yData)



    def alphaFunc(x, A, onset, tau):

        f = A * ((x - onset) / tau) * np.exp(-(x - onset - tau) / tau)
        f[x < onset] = 0

        return f

    def accept_test(f_new, x_new, f_old, x_old):
        [A1, onset1, tau1, A2, onset2, tau2, offset] = x_new

        return all([x > 0 for x in [A1, A2, tau1, tau2]]) and bool(onset1 >= -10) and bool(onset1 <= 15) \
               and bool(onset2 >= 0) and bool(onset2 <= 40) \
               and bool(offset <= 2) and bool(offset >= -2)

        # return all([x >= 0 for x in [A1, A2, tau1, tau2]])


    def take_step(x_old):


        rand = np.random.rand(*np.shape(x_old))
        # rand[1] = 2 * rand[1] - 1
        # rand[4] = 2 * rand[4] - 1
        # rand[6] = 2 * rand[6] - 1
        randT = 2 * rand - 1
        stepSizes = np.array([10, 5, 50, 10, 5, 50, 0.4])
        xNew = x_old + stepSizes * randT
        # temp = xNew[[0, 2, 3, 4, 5]]
        # temp[temp <= 0] = 0.5
        # xNew[[0, 2, 3, 5]] = temp
        # xNew[1] = min(15, max(xNew[1], -10))
        # xNew[4] = min(40, max(xNew[4], 0))
        # xNew[6] = min(2, max(xNew[6], -2))
        # assert accept_test(None, xNew, None, None), 'Step out of bounds'
        return xNew

    def diffOfAlphafuncs(x, A1, onset1, tau1, A2, onset2, tau2, offset):

        f1 = alphaFunc(x, A1, onset1, tau1)
        f2 = alphaFunc(x, A2, onset2, tau2)

        return f1 - f2 + offset

    def least_sq(pars):

        fitSignal = diffOfAlphafuncs(xData, *pars)

        return np.linalg.norm(yData - fitSignal)


    optRes = basinhopping(least_sq, x0=p0, T=10,
                            # niter_success=50,
                            niter=500,
                            # callback=callback,
                            take_step=take_step,
                            # accept_test=accept_test,
                            minimizer_kwargs={'method': 'L-BFGS-B'},
                            # minimizer_kwargs={'method': 'Powell'},
                            # disp=True
                            )

    print(optRes.x, optRes.fun)

    return optRes.x, optRes.fun


pool = Pool(processes=6)
argGen = ArgGenIterator(pars['xData'], pars['yData'], pars['p0s'])
vals = pool.map(fitAlpha, argGen)
xs = [x[0].tolist() for x in vals]
fs = [x[1].tolist() for x in vals]

outFile = pars['outFile']
with open(outFile, 'w') as fle:
    json.dump({'p0s': pars['p0s'], 'xs': xs, 'fs': fs}, fle)

