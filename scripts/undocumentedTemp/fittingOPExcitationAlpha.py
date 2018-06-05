from GJEMS.ephys.rawDataAnalyse import RawDataAnalyser
from scipy.signal import iirfilter, freqz, lfilter, kaiserord, firwin
import numpy as np
from neo import AnalogSignal
from matplotlib import pyplot as plt
plt.ion()
import os
import quantities as qu
import time
import json
import operator
import subprocess
from scipy.signal import medfilt

# enable_debugging = True
enable_debugging = False

if enable_debugging:
    import ipdb


mplPars = { 'text.usetex'       :    True,
            'axes.labelsize'    :   'large',
            'font.family'       :   'serif',
            'font.sans-serif'   :   'computer modern roman',
            }

for a, b in mplPars.items():
            plt.rcParams[a] = b


def LPFilterKaiser(signal, cutoff=100, transitionWidth=40, rippleDB=20):

    cutoff *= qu.Hz
    nyqFreq = signal.sampling_rate / 2


    transitionWidth = transitionWidth * qu.Hz

    N, beta = kaiserord(rippleDB, transitionWidth / nyqFreq)

    tapsLP = firwin(N, cutoff/nyqFreq, window=('kaiser', beta))

    delay = (N - 1) * 0.5 * signal.sampling_period


    filteredSignal = AnalogSignal(
                                           signal=lfilter(tapsLP, 1.0, signal.magnitude),
                                           sampling_rate=signal.sampling_rate,
                                           units=signal.units,
                                           t_start=signal.t_start - delay
                                          )


    return delay, filteredSignal


def HPFilterKaiser(signal, cutoff=10, transitionWidth=40, rippleDB=20):

    cutoff *= qu.Hz
    nyqFreq = signal.sampling_rate / 2

    transitionWidth = transitionWidth * qu.Hz

    N, beta = kaiserord(rippleDB, transitionWidth / nyqFreq)

    tapsLP = firwin(N, cutoff/nyqFreq, window=('kaiser', beta))

    temp = np.zeros((N,))
    temp[(N - 1) / 2] = 1
    tapsHP = temp - tapsLP

    delay = (N - 1) * 0.5 * signal.sampling_period


    filteredSignal = AnalogSignal(
                                           signal=lfilter(tapsHP, 1.0, signal.magnitude),
                                           sampling_rate=signal.sampling_rate,
                                           units=signal.units,
                                           t_start=signal.t_start - delay
                                          )


    return delay, filteredSignal


def getNoiseVar(resp):

    temp = resp - np.median(resp)
    return np.median(np.abs(temp)) / 0.6745


def getSpikeAmps(resp, spikeTimes):

    spikeInds = map(int, (spikeTimes - resp.t_start) * resp.sampling_rate)
    return resp[spikeInds]

def alphaFunc(x, A, onset, tau):

    f = A * ((x - onset) / tau) * np.exp(-(x - onset - tau) / tau)
    f[x < onset] = 0

    return f


def diffOfAlphafuncs(x, A1, onset1, tau1, A2, onset2, tau2):

    f1 = alphaFunc(x, A1, onset1, tau1)
    f2 = alphaFunc(x, A2, onset2, tau2)

    return f1 - f2





def getLocalMaxima(a):
    a = np.array(a)
    b = np.arange(a.shape[0])
    c = np.hstack((False, a[1:] > a[:-1])) & np.hstack((a[:-1] > a[1:], False))
    return b[c], a[c]




expNames = [
            '130313-4Rh',
            '130322-1LY',
            '130326-2Rh',
            '130408-1LY',
            '130425-1Al',
            '130501-2Rh',
            '130523-3LY',
            '130605-1LY',
            '130605-2LY',
            '130705-1LY',
            '140424-1LY',
            '140701-1Al',
            '140813-3Al',
            '140930-1Al',
            '140917-1Al',
            '141030-1Al',
            ]

homeFolder = os.path.expanduser('~')

dirpath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFittingAlpha/')
toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')
scalesFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/scales.json')
tempParFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting/tempPars265All.json')

if not os.path.exists(resDir):
    os.mkdir(resDir)
with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

with open(scalesFile, 'r') as fle:
    scales = json.load(fle)

Ts = 4.8e-5 * qu.s

traceLength = 1 * qu.s
# nPts = int(traceLength / Ts)

traceStart = 0 * qu.s

refInd = 0

medianFilterDur = 80 * qu.ms

nIterations = 200


types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [1 * qu.s, 1 * qu.s, 1 * qu.s]


for expName in expNames:

    print('Doing ' + expName)

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


    for resp in resps[265.0]:



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

            # tempNoiseVar = getNoiseVar(temp)
            # presSigScaled[presSigScaled > 5 * tempNoiseVar] = 5 * tempNoiseVar
            # presSigScaled[presSigScaled < -5 * tempNoiseVar] = -5 * tempNoiseVar

            # delay, presSigScaledFiltered = LPFilterKaiser(presSigScaled, cutoff=15, transitionWidth=10, rippleDB=60)
            # presSigScaledFiltered = medfilt(presSigScaled, medianFilterLen)

            respNormedSigs.append(presSigScaled)
            # respNormedSigsFiltered.append(presSigScaledFiltered)

        if len(respNormedSigs) == 3:

            # respFBFiltered = AnalogSignal(signal=np.concatenate(respNormedSigsFiltered),
            #                       units=qu.dimensionless,
            #                       sampling_rate=temp.sampling_rate,
            #                       t_start=-typeDurs[0],
            #                       )
            # respFB = np.concatenate(respNormedSigs)
            #
            # normedFilteredSigs.append(respFBFiltered)
            # normedSigs.append(respFB)

            respFBSignal = np.concatenate(respNormedSigs)
            baseLineSig = np.concatenate((respNormedSigs[0], respNormedSigs[2]))
            respFBSignal = respFBSignal - np.median(baseLineSig)
            # tempVarNoise = getNoiseVar(baseLineSig)
            # respFBSignal[respFBSignal > 5 * tempVarNoise] = 5 * tempVarNoise
            # respFBSignal[respFBSignal < -5 * tempVarNoise] = -5 * tempVarNoise

            respFB = AnalogSignal(signal=respFBSignal,
                                  units=qu.mV,
                                  sampling_rate=temp.sampling_rate,
                                  t_start=-typeDurs[0],
                                  )
            delay, respFBFiltered = LPFilterKaiser(respFB, cutoff=50, transitionWidth=40, rippleDB=60)

            # respFBFilteredSig = medfilt(respFBSignal, medianFilterLen)
            # respFBFiltered = AnalogSignal(signal=respFBFilteredSig,
            #                       units=qu.mV,
            #                       sampling_rate=temp.sampling_rate,
            #                       t_start=-typeDurs[0],
            #                       )
            normedFilteredSigs.append(respFBFiltered)
            normedSigs.append(respFB)

    avgSig = reduce(lambda x, y: (x + y), normedSigs) / len(normedSigs)
    avgSigLP = reduce(lambda x, y: (x + y), normedFilteredSigs) / len(normedFilteredSigs)
    # delay, avgSigLP = LPFilterKaiser(avgSig, cutoff=15, transitionWidth=10, rippleDB=60)


    signal2FitFull = avgSigLP
    fitStartT = (traceStart - signal2FitFull.t_start) * signal2FitFull.sampling_rate
    fitStart = int(fitStartT.simplified)
    fitEndT = (traceStart + traceLength - signal2FitFull.t_start) * signal2FitFull.sampling_rate
    fitEnd = int(fitEndT.simplified)

    signal2Fit = signal2FitFull[fitStart: fitEnd + 1]

    xDataT = signal2Fit.times
    xDataT.units = qu.ms
    xData = xDataT.magnitude
    yData = signal2Fit.magnitude

    def getp0s(nSamples, leftEnds, widths):

        assert len(leftEnds) == len(widths)
        all = np.empty((nSamples, len(leftEnds)))
        for ind in range(len(leftEnds)):
            all[:, ind] = np.random.rand(nSamples) * widths[ind] + leftEnds[ind]

        return all.tolist()

    # --------------------------------------------------------------------------------------------
    # A0, onset0, tau0, A1, onset1, tau1

    leftEnds = [10, 0, 50, 10, 0, 50, -1]
    widths = [50, 10, 200, 50, 10, 200, 2]

    p0s = getp0s(nIterations, leftEnds, widths)
    # print(p0s)

    outFile = os.path.join(resDir, expName + '.json')
    with open(tempParFile, 'w') as fle:
        json.dump({'p0s': p0s, 'outFile': outFile, 'xData': xData.tolist(), 'yData': yData.tolist(),
                   }, fle)

    subprocess.call(['python', 'scripts/fitAlphaFunc.py', tempParFile])
    # --------------------------------------------------------------------------------------------

    with open(outFile, 'r') as fle:
        pars = json.load(fle)
        p0s = pars['p0s']
        xs = pars['xs']
        fs = pars['fs']




















