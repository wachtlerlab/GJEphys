import os
import matplotlib.pyplot as plt
import quantities as qu
import json
from GJEMS.ephys.rawDataAnalyse import RawDataAnalyser
import operator
import numpy as np
from neo import AnalogSignal
from GJEMS.ephys.NEOFuncs import sliceAnalogSignal

plt.ion()
mplPars = {'text.usetex': True,
           'axes.labelsize': 'large',
           'font.family': 'sans-serif',
           'font.sans-serif': 'computer modern roman',
           'font.size': 24,
           'font.weight': 'black',
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           }

for x, y in mplPars.iteritems():
    plt.rcParams[x] = y

def accept_test1(x_new):
    [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = x_new

    return all(x > 0 for x in [taur1, taud1, taur2, taud2, A1, A2]) and bool(abs(offset) < 3) \
            and bool(-50 <= onset1 <= 50) and bool(0 <= onset2 <= 100) \
            and all(x < 100 for x in [taur1, taud1, taur2]) and bool(taud2 < 1000)

def doubleExpFun(x, Ar, Ad, t0, itaur, itaud):

        expd = Ad * np.exp(-itaud * (x - t0))
        expr = Ar * np.exp(-itaur * (x - t0))

        doubleExp = expd - expr
        doubleExp[x < t0] = (Ad - Ar)

        return doubleExp

expNames = [
    '130313-4Rh',
    '130322-1LY',
    '130326-2Rh',
    '130408-1LY',
    '130425-1Al',
    '130501-2Rh',
    '130705-1LY',
    '140424-1LY',

    '130523-3LY',
    '130605-1LY',
    '130605-2LY',
    '140701-1Al',
    '140813-3Al',
    '140930-1Al',
    '140917-1Al',
    '141030-1Al',
]

bestTrials = [16, 3, 2, None, 1, 4, 0, None, 0, None, 2, None, 0, None, 2, 2]
traceLength = 1.2 * qu.s

traceStart = -0.1 * qu.s

types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [1 * qu.s, 1 * qu.s, 1 * qu.s]

homeFolder = os.path.expanduser('~')

dirpath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting2Exp/')
toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')

fig, ax = plt.subplots(figsize=(14, 11.2))

with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

plottedExpNames = []
for expInd, expName in enumerate(expNames):

    print('Doing ' + expName)
    nrnResDir = os.path.join(resDir, expName)

    intervalToIgnore = None
    if expName in toIgnore:
        intervalToIgnore = toIgnore[expName]

    def shouldBeIgnored(resp):

        if intervalToIgnore is None:
            return False
        else:
            respInIgnoreIntervals = [(x * qu.s <= resp.t_start <= y * qu.s) | (x * qu.s <= resp.t_stop <= y * qu.s)
                                        for x, y in intervalToIgnore]
            return reduce(operator.or_, respInIgnoreIntervals)

    rda = RawDataAnalyser(expName, dirpath)

    resps = rda.getContResps(types=types, freqs=[265.0])

    normedFilteredSigs = []
    normedSigs = []
    t_starts = []
    nrn_fs = []

    for respInd, resp in enumerate(resps[265.0]):

        print('Doing Trial ' + str(respInd))

        t_starts.append(resp['DuringStimulus'].t_start)

        respNormedSigs = []
        respNormedSigsFiltered = []

        for typeInd, tpye in enumerate(types):

            temp = resp[tpye]
            if shouldBeIgnored(temp):
                break
            centeredSig = temp.magnitude
            # centeredSig = temp1 - np.median(temp1)
            typeLen = int(typeDurs[typeInd] * temp.sampling_rate)
            presSigScaled = np.zeros((typeLen, ))
            sigLen = min(typeLen, temp.shape[0])

            # medianFilterLen = int((medianFilterDur * temp.sampling_rate).simplified)
            # if not medianFilterLen % 2:
            #     medianFilterLen += 1

            if tpye == 'BeforeStimulus':

                presSigScaled[-sigLen:] = centeredSig[-sigLen:]
                if sigLen < typeLen:
                    presSigScaled[:-sigLen] = np.median(centeredSig[-sigLen:])

            elif tpye == 'DuringStimulus':

                if sigLen < typeLen:
                    break
                presSigScaled[:sigLen] = centeredSig[:sigLen]
            else:
                if sigLen < typeLen:
                    presSigScaled[sigLen:] = np.median(centeredSig[:sigLen])
                presSigScaled[:sigLen] = centeredSig[:sigLen]


            respNormedSigs.append(presSigScaled)


        if len(respNormedSigs) == 3:

            respFBSignal = np.concatenate(respNormedSigs)
            baseLineSig = np.concatenate((respNormedSigs[0], respNormedSigs[2]))
            respFBSignal = respFBSignal - np.median(baseLineSig)

            respFB = AnalogSignal(signal=respFBSignal,
                                  units=qu.mV,
                                  sampling_rate=temp.sampling_rate,
                                  t_start=-typeDurs[0],
                                  )

            normedSigs.append(respFB)

        signal2Fit = sliceAnalogSignal(respFB, traceStart, traceStart + traceLength)

        xDataT = signal2Fit.times
        xDataT.units = qu.ms
        xData = xDataT.magnitude

        outFile = os.path.join(nrnResDir, expName + '_Trial' + str(respInd) + '.json')
        if os.path.exists(outFile):

            with open(outFile, 'r') as fle:
                pars = json.load(fle)
                p0s = pars['p0s']
                xs = pars['xs']
                fs = pars['fs']

            xs_accepted = []
            fs_accepted = []
            acceptance = []

            for xInd, x in enumerate(xs):

                f = fs[xInd]

                [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset] = x

                if taud1 < taur1:
                    [A1, taud1, taur1] = [-A1, taur1, taud1]

                if taud2 < taur2:
                    [A2, taud2, taur2] = [-A2, taur2, taud2]

                if A1 < 0 and A2 < 0:
                    x = [-A2, -A1, t2, t1, taur2, taud2, taur1, taud1, offset]
                    xs[xInd] = x

                acc = accept_test1(x)
                acceptance.append(acc)

                if acc:
                    xs_accepted.append(x)
                    fs_accepted.append(f)

            if fs_accepted:

                bestFitInd = int(np.argmin(fs_accepted))
                xBest = xs_accepted[bestFitInd]
                [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
                d1 = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)
                d2 = -doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)
                exi2Inh = (-d1.max() / d2.min())
                print(exi2Inh)
                if 0.1 <= exi2Inh <= 1:
                    nrn_fs.append(min(fs_accepted))
                else:
                    nrn_fs.append(np.inf)

            else:
                nrn_fs.append(min(fs))

        elif respInd < len(normedSigs):
            normedSigs.pop(respInd)

    # if nrn_fs:
    #     bestRespInd = int(np.argmin(nrn_fs))
    #     plottedExpNames.append(expName)

    bestRespInd = bestTrials[expInd]
    if not bestRespInd is None:
        print(bestRespInd)
        plottedExpNames.append(expName)

        bestRespFull = normedSigs[bestRespInd]

        bestResp = sliceAnalogSignal(bestRespFull, traceStart, traceStart + traceLength)

        bestRespMin = bestResp.min()
        dynamicRange = bestResp.max() - bestRespMin

        bestResp2Plot = ((bestResp - bestRespMin) * 100 * bestResp.units / dynamicRange) + bestRespMin

        ax.plot(bestResp2Plot.times, bestResp2Plot.magnitude + 120 * (len(plottedExpNames) - 1))


ax.set_xlabel('Time (s)')
ax.set_yticks(np.arange(0, len(plottedExpNames)) * 120)
ax.set_yticklabels(plottedExpNames)
fig.tight_layout()
fig.canvas.draw()




