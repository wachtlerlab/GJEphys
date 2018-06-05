import json
import quantities as qu
import numpy as np
import os
from matplotlib import pyplot as plt
import seaborn as sns
plt.ion()
homeFolder = os.path.expanduser('~')

mplPars = {'text.usetex': True,
           'axes.labelsize': 'large',
           'font.family': 'sans-serif',
           'font.sans-serif': 'computer modern roman',
           'font.size': 24,
           'font.weight': 'black',
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           }

sns.set(rc=mplPars)

def accept_test1(x_new):
    [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = x_new

    return all(x > 0 for x in [taur1, taud1, taur2, taud2, A1, A2]) and bool(abs(offset) < 3) \
            and bool(-50 <= onset1 <= 50) and bool(0 <= onset2 <= 100) \
            and all(x < 100 for x in [taur1, taud1, taur2])

def doubleExpFun(xSig, Ar, Ad, t0, itaur, itaud):

    expd = Ad * np.exp(-itaud * (xSig - t0))
    expr = Ar * np.exp(-itaur * (xSig - t0))

    doubleExp = expd - expr
    doubleExp[xSig < t0] = (Ad - Ar)

    return doubleExp

def twoDoubleExps(xSig, A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset):

        d1 = doubleExpFun(xSig, A1, A1, t1, 1 / taur1, 1 / taud1)
        d2 = doubleExpFun(xSig, A2, A2, t2, 1 / taur2, 1 / taud2)

        return d1 - d2 + offset

Ts = 4.8e-5 * qu.s
Ts.units = qu.ms
sigDur = 1 * qu.s
sigDur.units = qu.ms

xData = np.arange(0, sigDur.magnitude, Ts.magnitude)

expNames = [
    '130313-4Rh',
    '130322-1LY',
    '130326-2Rh',
    '130408-1LY',
    '130425-1Al',
    '130501-2Rh',
    '130705-1LY',
    # '140424-1LY',
    '130523-3LY',
    '130605-1LY',
    '130605-2LY',
    # '140701-1Al',
    '140813-3Al',
    '140930-1Al',
    '140917-1Al',
    '141030-1Al',
]

# doeParsDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Backups/OPExcitationFitting_2exp_1/')
doeParsDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting2Exp/')
outJSONFile = os.path.join(doeParsDir, 'bestPars.json')

Ts = 4.8e-5 * qu.s
Ts.units = qu.ms
sigDur = 1 * qu.s
sigDur.units = qu.ms

xData = np.arange(0, sigDur.magnitude, Ts.magnitude).tolist()

info = "These parameters can be used to generate the assumed input current waveforms for DL-Int-1 neurons. " \
       "Please have a look at modelling.pdf and vizFittingOPExcitation2Exp.py, especially the function" \
       "twoDoubleExps() and the definition of fitSignal on line 455 for more info"

outJSON = {'xData': xData, 'bestPars': {}, 'info': info}

for expInd, expName in enumerate(expNames):

    print('Doing ' + expName)

    outFile = os.path.join(doeParsDir, expName + '.json')
    with open(outFile, 'r') as fle:
        pars = json.load(fle)
        p0s = pars['p0s']
        xs = pars['xs']
        fs = pars['fs']

    xs_accepted = []
    fs_accepted = []
    acceptance = []

    for x, f in zip(xs, fs):

        [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset] = x

        if taud1 < taur1:
            [A1, taud1, taur1] = [-A1, taur1, taud1]

        if taud2 < taur2:
            [A2, taud2, taur2] = [-A2, taur2, taud2]

        if A1 < 0 and A2 < 0:
            x = [-A2, -A1, t2, t1, taur2, taud2, taur1, taud1, offset]

        acc = accept_test1(x)
        acceptance.append(acc)

        if acc:
            xs_accepted.append(x)
            fs_accepted.append(f)

    if xs_accepted:
        xBest = xs_accepted[int(np.argmin(fs_accepted))]
        outJSON['bestPars'][expName] = xBest


with open(outJSONFile, 'w') as fle:
    json.dump(outJSON, fle)

