from GJEMS.ephys.rawDataAnalyse import RawDataAnalyser
from scipy.signal import iirfilter, freqz, lfilter, kaiserord, firwin
import numpy as np
from neo import AnalogSignal
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize, basinhopping
plt.ion()
import os
import quantities as qu
import time
import json
import operator
import subprocess

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



def doubleExpFun(x, Ar, Ad, t0, itaur, itaud):

    expd = Ad * np.exp(-itaud * (x - t0))
    expr = Ar * np.exp(-itaur * (x - t0))

    doubleExp = expd - expr
    doubleExp[x < t0] = (Ad - Ar)

    return doubleExp





def threeDoubleExps(x, Ar1, Ad1, Ar2, Ad2, Ar3, Ad3, t1, t2, t3, itaur1, itaud1, itaur2, itaud2, itaur3, itaud3):

    if np.any(np.array([Ar1, Ad1, Ar2, Ad2, Ar3, Ad3, itaur1, itaud1, itaur2, itaud2, itaur3, itaud3]) < 0):

        return np.ones_like(x) * 100
    #
    # elif (Ad1 < Ar1) or (Ad2 > Ar2) or (Ad3 > Ar3):
    #     return np.ones_like(x) * 100

    return doubleExpFun(x, Ar1, Ad1, t1, itaur1, itaud1) \
               + doubleExpFun(x, Ar2, Ad2, t2, itaur2, itaud2) \
               + doubleExpFun(x, Ar3, Ad3, t3, itaur3, itaud3)

def accept_test(f_new, x_new, f_old, x_old):

    [Ar1, Ad1, Ar2, Ad2, Ar3, Ad3, t1, t2, t3, itaur1, itaud1, itaur2, itaud2, itaur3, itaud3] = x_new

    isPos = lambda x: x > 0

    logicalAnd = lambda x, y: x & y

    areAllPos = reduce(logicalAnd,
                  map(isPos, [Ar1, Ad1, Ar2, Ad2, Ar3, Ad3, itaur1, itaud1, itaur2, itaud2, itaur3, itaud3]))

    areTimesAcceptable = reduce(logicalAnd, map(lambda x: -5 < x < 10, [t1, t2, t3]))



    decision = bool(areAllPos
                    and (itaud2 < itaur2) and (itaud3 < itaur3)
                    and (t1 <= t2) and (t2 <= t3)
                    and (itaud1 > 0.3) and (itaud2 > 0.3) and (itaur1 > 0.3)
                    )
    if decision:
        print(f_new, decision)
    return decision


def getLocalMaxima(a):
    a = np.array(a)
    b = np.arange(a.shape[0])
    c = np.hstack((False, a[1:] > a[:-1])) & np.hstack((a[:-1] > a[1:], False))
    return b[c], a[c]




# expName = '130705-1LY'
expName = '130425-1Al'



dirpath = '/home/ajay/DataAndResults/GJEphys/NIXFiles/'
resDir = '/home/ajay/DataAndResults/GJEphys/OPExcitationFitting/'
toIgnoreFile = '/home/ajay/DataAndResults/GJEphys/NIXFiles/toIgnore.json'
tempParFile = '/home/ajay/DataAndResults/GJEphys/OPExcitationFitting/tempParsFreq.json'
with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

Ts = 4.8e-5 * qu.s

traceLength = 1
nPts = int(traceLength / Ts)

traceStart = 0 * qu.s

# freqs = [100]
# freqs = [200]
# freqs = [265]
# freqs = [300]
# freqs = [400]
freqs = [100, 200, 265, 300, 400]


types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [3 * qu.s, 1 * qu.s, 3 * qu.s]

print('Doing ' + expName)
for freq in freqs:
    print('Doing ' + str(freq) + 'Hz')


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

    resps = rda.getContResps(types=types, freqs=freqs)

    normedFilteredSigs = []

    refSig = resps[freq][0]['DuringStimulus']
    refSigCentered = refSig - np.median(refSig)

    for resp in resps[freq]:



        respNormedSigs = []

        for typeInd, tpye in enumerate(types):

            temp = resp[tpye]
            if shouldBeIgnored(temp):
                break
            temp1 = temp.magnitude
            centeredSig = temp1 - np.median(temp1)
            typeLen = int(typeDurs[typeInd] * temp.sampling_rate)
            presSigScaled = np.zeros((typeLen, ))
            sigLen = min(typeLen, temp.shape[0])
            tempNoiseVar = getNoiseVar(temp)

            if tpye == 'BeforeStimulus':

                presSigScaled[-sigLen:] = centeredSig[-sigLen:]

            else:
                presSigScaled[:sigLen] = centeredSig[:sigLen]

            presSigScaled[presSigScaled > 50 * tempNoiseVar] = 50 * tempNoiseVar
            presSigScaled[presSigScaled < -50 * tempNoiseVar] = -50 * tempNoiseVar

            respNormedSigs += [presSigScaled]

        if len(respNormedSigs) == 3:

            respFB = AnalogSignal(signal=np.concatenate(respNormedSigs),
                                  units=qu.dimensionless,
                                  sampling_rate=temp.sampling_rate,
                                  t_start=-3 * qu.s,
                                  )
            normedFilteredSigs.append(respFB)


    avgSig = reduce(lambda x, y: (x + y), normedFilteredSigs) / len(normedFilteredSigs)
    delay, avgSigLP = LPFilterKaiser(avgSig, cutoff=15, transitionWidth=10, rippleDB=60)

    # fig, ax = plt.subplots(figsize=(10, 8))
    # ax.plot(avgSig.times, avgSig, 'b')
    # ax.plot(avgSigLP.times, avgSigLP, 'r')
    # ax.legend(['Sig', 'SigLP'])


    delayms = delay.copy()
    delayms.units = qu.ms
    delayms = delayms.magnitude
    print(delayms)

    def twoDoubleExps(x, A1, A2, t1, t2, taur1, taud1, taur2, taud2):

        d1 = doubleExpFun(x, A1, A1, t1, 1 / taur1, 1 / taud1)
        d2 = doubleExpFun(x, A2, A2, t2, 1 / taur2, 1 / taud2)

        # if np.any(np.array([A1, A2, t1, t2, itaur1, itaud1, itaur2, itaud2]) < 0)\
        #         or np.any(np.array([t1, t2]) < delayms - 10):
        #
        #     return np.ones_like(x) * 100
        # to avoid inverted double exponential for the lower delay one
        # elif (Ad1 < Ar1):
        #     return np.zeros_like(x)
        # else:
        return d1 - d2


    signal2FitFull = avgSigLP
    sigStart = traceStart - delay
    totalPoints = nPts + int(delay * signal2FitFull.sampling_rate)
    t = sigStart + np.arange(totalPoints) * signal2FitFull.sampling_period
    fitStart = int((sigStart - signal2FitFull.t_start) * signal2FitFull.sampling_rate)

    # fig, ax = plt.subplots(figsize=(10, 8))
    # ax.plot(avgSig.times, avgSig, 'b')
    # ax.plot(avgSigLP.times, avgSigLP, 'r')
    # ax.legend(['Sig', 'SigLP'])
    # ax.set_title(str(freq))
    #
    # fig.canvas.draw()


    baseLine = np.concatenate((signal2FitFull[0:fitStart], signal2FitFull[fitStart + totalPoints + 1:]), axis=1).mean()
    signal2Fit = signal2FitFull[fitStart: fitStart + totalPoints + 1] - baseLine
    print(baseLine)
    xData = (signal2Fit.times - signal2Fit.t_start).magnitude * 1e3
    yData = signal2Fit.magnitude

    f = lambda x: np.concatenate((np.random.rand(2) * 10,
                                  np.random.rand(2) * 100 + delayms - 50,
                                  np.random.rand(4) * 50 + 50),
                                 axis=1).tolist()
    p0s = map(f, range(10))
    outFile = os.path.join(resDir, expName + str(int(freq)) + '.npz')

    with open(tempParFile, 'w') as fle:
        json.dump({'p0s': p0s, 'outFile': outFile, 'xData': xData.tolist(), 'yData': yData.tolist(),
                   'filterDelay': float(delayms)}, fle)

    subprocess.call(['python', 'scripts/fitDoubleExp.py', tempParFile])

    with open(outFile, 'r') as fle:
        pOpt = json.load(fle)

    [A1, A2, t1, t2, taur1, taud1, taur2, taud2] = pOpt
    print(np.round([A1, A2], 4).T)
    print(np.round([t1, t2], 4).T)
    print(np.round([taur1, taud1, taur2, taud2], 4).T)
    fitSignal = twoDoubleExps(xData, *pOpt)

    if t2 > t1:
        doubleExp1 = doubleExpFun(xData, A1, A1, t1, 1 / taur1, 1 / taud1)
        doubleExp2 = -doubleExpFun(xData, A2, A2, t2, 1 / taur2, 1 / taud2)


    else:
        doubleExp2 = doubleExpFun(xData, A1, A1, t1, 1 / taur1, 1 / taud1)
        doubleExp1 = -doubleExpFun(xData, A2, A2, t2, 1 / taur2, 1 / taud2)


    print(fitSignal.max(), fitSignal.min())
    print(doubleExp1.max(), doubleExp2.min())

    # fig1, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(10, 16))
    # # ax1.plot((centeredSignal.times - centeredSignal.t_start) * 1e3, centeredSignal, 'r', label='raw')
    # ax1[0].plot((signal2Fit.times - signal2Fit.t_start) * 1e3, signal2Fit, 'g', label='filtered')
    # ax1[0].plot(xData, fitSignal, 'b', label='fit')
    #
    # ax1[0].legend()
    #
    # ax1[1].plot(xData, fitSignal, 'b', label='fit')
    # ax1[1].plot(xData, doubleExp2, 'r', label='doubleExp2')
    # ax1[1].plot(xData, doubleExp1, 'c', label='doubleExp1')
    # ax1[1].legend()
    # ax1[0].set_title(str(freq))



    # fig1.canvas.draw()
    # raw_input()














