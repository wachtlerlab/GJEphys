from GJEphys.rawDataAnalyse import RawDataAnalyser
from scipy.signal import iirfilter, freqz, lfilter, kaiserord, firwin
import numpy as np
from neo import AnalogSignal
from matplotlib import pyplot as plt
from scipy.signal import medfilt
plt.ion()
import os
import quantities as qu
import time
import json
import operator

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


def getFirstZeroCrossingIndex(sig):

    sigN = np.array(sig)

    zc = np.where(np.logical_and(sigN[:-1] <= 0, sigN[1:] >= 0))[0]

    if len(zc):
        return zc[0]
    else:
        return None


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



def doubleExpFun(xSig, Ar, Ad, t0, itaur, itaud):

    expd = Ad * np.exp(-itaud * (xSig - t0))
    expr = Ar * np.exp(-itaur * (xSig - t0))

    doubleExp = expd - expr
    doubleExp[xSig < t0] = (Ad - Ar)

    return doubleExp

def twoDoubleExps(xSig, A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset):

        d1 = doubleExpFun(xSig, A1, A1, t1, 1 / taur1, 1 / taud1)
        d2 = doubleExpFun(xSig, A2, A2, t2, 1 / taur2, 1 / taud2)

        # if np.any(np.array([A1, A2, t1, t2, itaur1, itaud1, itaur2, itaud2]) < 0)\
        #         or np.any(np.array([t1, t2]) < delayms - 10):
        #
        #     return np.ones_like(x) * 100
        # to avoid inverted double exponential for the lower delay one
        # elif (Ad1 < Ar1):
        #     return np.zeros_like(x)
        # else:
        return d1 - d2 + offset


def alphaFunc(x, A, onset, tau):

    f = A * ((x - onset) / tau) * np.exp(-(x - onset - tau) / tau)
    f[x < onset] = 0

    return f


def diffOfAlphafuncs(x, A1, onset1, tau1, A2, onset2, tau2, offset):

    f1 = alphaFunc(x, A1, onset1, tau1)
    f2 = alphaFunc(x, A2, onset2, tau2)

    return f1 - f2 + offset


def threeDoubleExps(x, Ar1, Ad1, Ar2, Ad2, Ar3, Ad3, t1, t2, t3, itaur1, itaud1, itaur2, itaud2, itaur3, itaud3):

    if np.any(np.array([Ar1, Ad1, Ar2, Ad2, Ar3, Ad3, itaur1, itaud1, itaur2, itaud2, itaur3, itaud3]) < 0):

        return np.ones_like(x) * 100
    #
    # elif (Ad1 < Ar1) or (Ad2 > Ar2) or (Ad3 > Ar3):
    #     return np.ones_like(x) * 100

    return doubleExpFun(x, Ar1, Ad1, t1, itaur1, itaud1) \
               + doubleExpFun(x, Ar2, Ad2, t2, itaur2, itaud2) \
               + doubleExpFun(x, Ar3, Ad3, t3, itaur3, itaud3)



def getLocalMaxima(a):
    a = np.array(a)
    b = np.arange(a.shape[0])
    c = np.hstack((False, a[1:] > a[:-1])) & np.hstack((a[:-1] > a[1:], False))
    return b[c], a[c]


def accept_test1(x_new):
    [A1, onset1, tau1, A2, onset2, tau2, offset] = x_new

    return all(x > 0 for x in [tau1, tau2, A1, A2]) and bool(abs(offset) < 3) \
            and bool(-50 <= onset1 <= 50) and bool(0 <= onset2 <= 100)




expNames = [
            # '130313-4Rh',
            # '130322-1LY',
            # '130326-2Rh',
            # '130408-1LY',
            # '130425-1Al',
            # '130501-2Rh',
            # '130523-3LY',
            # '130605-1LY',
            # '130605-2LY',
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
with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

with open(scalesFile, 'r') as fle:
    scales = json.load(fle)

Ts = 4.8e-5 * qu.s

traceLength = 1 * qu.s
# nPts = int(traceLength / Ts)

traceStart = 0 * qu.ms

refInd = 0

medianFilterDur = 80 * qu.ms


types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [1 * qu.s, 1 * qu.s, 1 * qu.s]

fig0, ax0 = plt.subplots(figsize=(10, 8))
fig, ax = plt.subplots(figsize=(10, 8))
fig1, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(10, 8))

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
            # tempMed = np.median(centeredSig)
            # centeredSig = temp1 - tempMed
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



    ax.clear()
    avgSig = reduce(lambda x, y: (x + y), normedSigs) / len(normedSigs)
    avgSigLP = reduce(lambda x, y: (x + y), normedFilteredSigs) / len(normedFilteredSigs)

    avgSigVar = getNoiseVar(avgSig.magnitude) * avgSig.units
    avgSig[avgSig < -20 * avgSigVar] = -20 * avgSigVar
    avgSig[avgSig > 20 * avgSigVar] = 20 * avgSigVar

    # delay, avgSigLP = LPFilterKaiser(avgSig, cutoff=15, transitionWidth=10, rippleDB=60)
    ax.plot(avgSig.times, avgSig, 'b', label='Sig')
    ax.plot(avgSigLP.times, avgSigLP, 'r', label='SigLP')
    ax.set_ylabel('Voltage (mV)')
    ax.set_xlabel('time (s)')


    dynamicRange = min(np.max(normedSigs) - np.min(normedSigs), 25 * avgSigVar.magnitude)

    ax0.clear()
    for sigInd, sig in enumerate(normedFilteredSigs):
        temp = normedSigs[sigInd]
        temp[temp < -20 * avgSigVar] = -20 * avgSigVar
        temp[temp > 20 * avgSigVar] = 20 * avgSigVar
        tempM = temp.magnitude
        sigM = sig.magnitude
        ax0.plot(temp.times, tempM - np.median(tempM) - sigInd * dynamicRange, 'r')
        ax0.plot(sig.times, sigM - np.median(sigM) - sigInd * dynamicRange, 'b')

    fig0.canvas.draw()

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

    outFile = os.path.join(resDir, expName + '.json')
    with open(outFile, 'r') as fle:
        pars = json.load(fle)
        p0s = pars['p0s']
        xs = pars['xs']
        fs = pars['fs']

    xs_accepted = []
    fs_accepted = []
    acceptance = []

    for x, f in zip(xs, fs):

        [A1, onset1, tau1, A2, onset2, tau2, offset] = x

        if A1 < 0 and A2 < 0:
            x = [-A2, onset2, tau2, -A1, onset1, tau1, offset]

        acc = accept_test1(x)
        acceptance.append(acc)

        if acc:
            xs_accepted.append(x)
            fs_accepted.append(f)

    map(lambda x: x.clear(), ax1.flat)

    fig1.suptitle(expName + '(n=' + str(len(normedSigs)) + ')')
    if not xs_accepted:
        print('None of the solutions are acceptable. Plotting the best among the solutions')
        ax1[0].set_title('None of the solutions are acceptable. Plotting the best among the solutions')
        xBest = xs[int(np.argsort(fs)[0])]
    else:
        xBest = xs_accepted[int(np.argsort(fs_accepted)[0])]



    [A1, onset1, tau1, A2, onset2, tau2, offset] = xBest


    fitSignal = diffOfAlphafuncs(xData, *xBest)
    a1 = alphaFunc(xData, A1, onset1, tau1)
    a2 = -alphaFunc(xData, A2, onset2, tau2)
    fitError = np.linalg.norm(fitSignal - yData)
    #
    # print(fitSignal.max(), fitSignal.min())
    # print(alpha1.max(), alpha2.min())

    print(round(offset, 4))
    print(np.round([A1, onset1, tau1], 4).T)
    print(np.round([A2, onset2, tau2], 4).T)
    print(fitError)

    temp = float((qu.ms / signal2Fit.t_start.units).simplified)
    ax.plot(xData * temp, fitSignal, 'g', label='fit')
    ax.legend()
    fig.canvas.draw()

    # ax1.plot((centeredSignal.times - centeredSignal.t_start) * 1e3, centeredSignal, 'r', label='raw')
    ax1[0].plot(xData, yData, 'g', label='filtered')
    ax1[0].plot(xData, fitSignal, 'b', label='fit')
    ax1[0].text(xData[int(len(xData) * 0.15)], 0.75 * max(yData), 'SSE=' + str(round(fitError, 2)), fontsize=14)

    ax1[0].legend()
    ax1[0].grid(True)

    ax1[1].plot(xData, fitSignal, 'b', label='fit')
    ax1[1].plot(xData, a1, 'r', label='alpha1')
    ax1[1].plot(xData, a2, 'c', label='alpha2')
    ax1[1].grid(True)
    ax1[1].legend()
    ax1[1].set_ylabel('Voltage (mV)')
    ax1[1].set_xlabel('time (' + signal2Fit.t_start.units.dimensionality.string + ')')




    fig1.canvas.draw()
    fig1.savefig(os.path.join(resDir, expName + '.eps'), dpi=600)

    if enable_debugging:
        ipdb.set_trace()
    else:
        raw_input()
    # plt.close('all')














