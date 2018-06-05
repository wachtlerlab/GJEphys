import numpy as np
from scipy.optimize import least_squares
import sys
import json
from multiprocessing import Pool

assert len(sys.argv) == 2, 'Only one argument, the path of the swcfile expected, ' + str(len(sys.argv)) + 'found'

parFile = sys.argv[1]
with open(parFile, 'r') as fle:
    pars = json.load(fle)


class ArgGenIterator:

    def __init__(self, xData, yData, p0s, bounds):


        self.xData = xData
        self.yData = yData
        self.p0s = p0s
        self.bounds = bounds
        self.pointsDone = 0

    def __iter__(self):

        self.pointsDone = 0
        return self

    def next(self):

        if self.pointsDone < len(self.p0s):
            toReturn = (self.xData, self.yData, self.p0s[self.pointsDone], self.bounds)
            self.pointsDone += 1
            return toReturn
        else:
            raise StopIteration

def fitDoubleExp(y):

    xData, yData, p0, bounds = y
    xData = np.array(xData)
    yData = np.array(yData)

    from GJEphys.doubleExpFitting import twoDoubleExps

    def least_sq(pars, xData, yData):

        fitSignal = twoDoubleExps(xData, *pars)
        return np.linalg.norm(yData - fitSignal)

    optRes = least_squares(least_sq, x0=p0, args=(xData, yData), verbose=0, bounds=bounds)

    # print(optRes.x, optRes.fun[0])
    return optRes.x, optRes.fun[0]


pool = Pool(processes=6)
argGen = ArgGenIterator(pars['xData'], pars['yData'], pars['p0s'], pars['bounds'])
vals = pool.map(fitDoubleExp, argGen)
xs = [x[0].tolist() for x in vals]
fs = [x[1].tolist() for x in vals]

outFile = pars['outFile']
with open(outFile, 'w') as fle:
    json.dump({'p0s': pars['p0s'], 'xs': xs, 'fs': fs}, fle)

