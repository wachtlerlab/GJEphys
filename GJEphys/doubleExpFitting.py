import numpy as np
import subprocess
import json
import os

def doubleExpFun(x, Ar, Ad, t0, itaur, itaud):

        expd = Ad * np.exp(-itaud * (x - t0))
        expr = Ar * np.exp(-itaur * (x - t0))

        doubleExp = expd - expr
        doubleExp[x < t0] = (Ad - Ar)

        return doubleExp



def twoDoubleExps(x, A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset):

        d1 = doubleExpFun(x, A1, A1, t1, 1 / taur1, 1 / taud1)
        d2 = doubleExpFun(x, A2, A2, t2, 1 / taur2, 1 / taud2)

        return d1 - d2 + offset


def accept_test1(x_new):
    [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = x_new

    return all(x > 0 for x in [taur1, taud1, taur2, taud2, A1, A2])\
            and bool(-50 <= onset1 <= 20) and bool(0 <= onset2 <= 40) and bool(abs(offset) <= 3)\
            and all(x < 100 for x in [taur1, taud1, taur2]) and bool(taud2 < 1000)


def getp0s(nSamples, leftEnds, widths):

    assert len(leftEnds) == len(widths)
    all = np.empty((nSamples, len(leftEnds)))
    for ind in range(len(leftEnds)):
        all[:, ind] = np.random.rand(nSamples) * widths[ind] + leftEnds[ind]

    return all.tolist()


def runManyDoubleExpFittings(xData, yData, nIterations, outFile, tempParFile):
    leftEnds = [10, 10, -50, 0, 0, 0, 0, 0, -1]
    widths = [50, 50, 50, 100, 100, 100, 100, 1000, 2]

    bounds = tuple(zip([0, np.inf], [0, np.inf], [-50, 50], [0, 100],
                       [0, 100], [0, 100], [0, 100], [0, 1000], [-3, 3]))

    p0s = getp0s(nIterations, leftEnds, widths)
    # print(p0s)


    with open(tempParFile, 'w') as fle:
        json.dump({'p0s': p0s, 'outFile': outFile, 'xData': xData.tolist(), 'yData': yData.tolist(),
                   'bounds': bounds}, fle)

    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'fitDoubleExp.py'), tempParFile])
