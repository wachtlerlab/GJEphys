import os
from matplotlib import pyplot as plt
import json
import quantities as qu
import operator
from GJEphys.rawDataAnalyse import RawDataAnalyser
import numpy as np
from neo import AnalogSignal
from scipy.signal import kaiserord, firwin, lfilter
plt.ion()

homeFolder = os.path.expanduser('~')

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

def getNoiseVar(resp):

    temp = resp - np.median(resp)
    return np.median(np.abs(temp)) / 0.6745

def alphaFunc(x, A, onset, tau):

        f = A * ((x - onset) / tau) * np.exp(-(x - onset - tau) / tau)
        f[x < onset] = 0

        return f

def diffOfAlphafuncs(x, A1, onset1, tau1, A2, onset2, tau2, offset):

        f1 = alphaFunc(x, A1, onset1, tau1)
        f2 = alphaFunc(x, A2, onset2, tau2)

        return f1 - f2 + offset

def doubleExpFun(x, Ar, Ad, t0, itaur, itaud):

        expd = Ad * np.exp(-itaud * (x - t0))
        expr = Ar * np.exp(-itaur * (x - t0))

        doubleExp = expd - expr
        doubleExp[x < t0] = (Ad - Ar)

        return doubleExp

def twoDoubleExps(x, A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset):

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
            return d1 - d2 + offset

def accept_test_alpha(x_new):
    [A1, onset1, tau1, A2, onset2, tau2, offset] = x_new

    return all(x > 0 for x in [tau1, tau2, A1, A2]) and bool(abs(offset) < 3) \
            and bool(-50 <= onset1 <= 50) and bool(0 <= onset2 <= 100) \
            and all(x < 100 for x in [tau1])


def accept_test_doe(x_new):
    [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = x_new

    return all(x > 0 for x in [taur1, taud1, taur2, taud2, A1, A2]) and bool(abs(offset) < 3) \
            and bool(-50 <= onset1 <= 50) and bool(0 <= onset2 <= 100) \
            and all(x < 100 for x in [taur1, taud1, taur2])


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


alphaParsDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Backups/OPExcitationFitting_alpha_2/')
doeParsDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Backups/OPExcitationFitting_2exp_1/')

dirpath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Backups/' \
                      + os.path.split(os.path.split(alphaParsDir)[0])[1] + '_' +
                      os.path.split(os.path.split(doeParsDir)[0])[1])
toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')
scalesFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/scales.json')
tempParFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting/tempPars265All.json')

if not os.path.isdir(resDir):
    os.mkdir(resDir)

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
fig1, ax1 = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))

aics_alpha = []
bics_alpha = []

aics_doe = []
bics_doe = []

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

    outFile = os.path.join(alphaParsDir, expName + '.json')
    with open(outFile, 'r') as fle:
        pars = json.load(fle)
        p0s = pars['p0s']
        xs = pars['xs']
        fs = pars['fs']

    xs_accepted_alpha = []
    fs_accepted_alpha = []
    acceptance_alpha = []

    for x, f in zip(xs, fs):

        [A1, onset1, tau1, A2, onset2, tau2, offset] = x

        if A1 < 0 and A2 < 0:
            x = [-A2, onset2, tau2, -A1, onset1, tau1, offset]

        acc = accept_test_alpha(x)
        acceptance_alpha.append(acc)

        if acc:
            xs_accepted_alpha.append(x)
            fs_accepted_alpha.append(f)

    map(lambda x: x.clear(), ax1.flat)

    fig1.suptitle(expName + '(n=' + str(len(normedSigs)) + ')')
    if not xs_accepted_alpha:
        ax1[0, 0].set_title('None of the solutions are acceptable')
        xBest_alpha = float('nan')
        aic = np.inf
        bic = np.inf

    else:
        xBest_alpha = xs_accepted_alpha[int(np.argsort(fs_accepted_alpha)[0])]

        [A1, onset1, tau1, A2, onset2, tau2, offset] = xBest_alpha

        fitSignal = diffOfAlphafuncs(xData, *xBest_alpha)
        a1 = alphaFunc(xData, A1, onset1, tau1)
        a2 = -alphaFunc(xData, A2, onset2, tau2)
        fitError = np.linalg.norm(fitSignal - yData)
        #
        # print(fitSignal.max(), fitSignal.min())
        # print(alpha1.max(), alpha2.min())

        # print(round(offset, 4))
        # print(np.round([A1, onset1, tau1], 4).T)
        # print(np.round([A2, onset2, tau2], 4).T)
        # print(fitError)


        # ax1.plot((centeredSignal.times - centeredSignal.t_start) * 1e3, centeredSignal, 'r', label='raw')
        ax1[0, 0].plot(xData, yData, 'g', label='filtered')
        ax1[0, 0].plot(xData, fitSignal, 'b', label='fit')
        ax1[0, 0].text(xData[int(len(xData) * 0.15)], 0.75 * max(yData), 'SSE=' + str(round(fitError, 2)), fontsize=14)
        ax1[0, 0].set_title('A, onset, tau')
        ax1[0, 0].legend()
        ax1[0, 0].grid(True)

        ax1[1, 0].plot(xData, fitSignal, 'b', label='fit')
        ax1[1, 0].plot(xData, a1, 'r', label='alpha1')
        ax1[1, 0].plot(xData, a2, 'c', label='alpha2')
        ax1[1, 0].grid(True)
        ax1[1, 0].legend()
        ax1[1, 0].set_title(str([round(x, 2) for x in [A1, onset1, tau1]])
                            + '\n'
                            + str([round(x, 2) for x in [A2, onset2, tau2]]))
        ax1[1, 0].set_ylabel('Voltage (mV)')
        ax1[1, 0].set_xlabel('time (' + signal2Fit.t_start.units.dimensionality.string + ')')

        N = xData.shape[0]

        # https://www.researchgate.net/post/What_is_the_AIC_formula
        aic = N * np.log((fitError ** 2) / N) + 2 * 7

        bic = N * np.log((fitError ** 2) / N) + np.log(N) * 7

    aics_alpha.append(aic)
    bics_alpha.append(bic)

    outFile = os.path.join(doeParsDir, expName + '.json')
    with open(outFile, 'r') as fle:
        pars = json.load(fle)
        p0s = pars['p0s']
        xs = pars['xs']
        fs = pars['fs']

    xs_accepted_doe = []
    fs_accepted_doe = []
    acceptance_doe = []

    for xInd, x in enumerate(xs):

        f = fs[xInd]

        [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset] = x

        if taud1 < taur1:
            [A1, taud1, taur1] = [-A1, taur1, taud1]

        if taud2 < taur2:
            [A2, taud2, taur2] = [-A2, taur2, taud2]

        if A1 < 0 and A2 < 0:
            x = [-A2, -A1, t2, t1, taur2, taud2, taur1, taud1, offset]
        else:
            x = [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset]

        xs[xInd] = x

        acc = accept_test_doe(x)
        acceptance_doe.append(acc)

        if acc:
            xs_accepted_doe.append(x)
            fs_accepted_doe.append(f)


    if not xs_accepted_doe:
        ax1[0, 0].set_title('None of the solutions are acceptable')
        xBest_doe = float('nan')
        aic = np.inf
        bic = np.inf

    else:
        xBest_doe = xs_accepted_doe[int(np.argsort(fs_accepted_doe)[0])]

        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest_doe

        fitSignal = twoDoubleExps(xData, *xBest_doe)
        d1 = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)
        d2 = -doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)
        fitError = np.linalg.norm(fitSignal - yData)
        #
        # print(fitSignal.max(), fitSignal.min())
        # print(alpha1.max(), alpha2.min())

        # print(round(offset, 4))
        # print(np.round([A1, onset1, tau1], 4).T)
        # print(np.round([A2, onset2, tau2], 4).T)
        # print(fitError)

        # ax1.plot((centeredSignal.times - centeredSignal.t_start) * 1e3, centeredSignal, 'r', label='raw')
        ax1[0, 1].plot(xData, yData, 'g', label='filtered')
        ax1[0, 1].plot(xData, fitSignal, 'b', label='fit')
        ax1[0, 1].text(xData[int(len(xData) * 0.15)], 0.75 * max(yData), 'SSE=' + str(round(fitError, 2)),
                       fontsize=14)
        ax1[0, 1].set_title('A, onset, taur, taud')
        ax1[0, 1].legend()
        ax1[0, 1].grid(True)

        ax1[1, 1].plot(xData, fitSignal, 'b', label='fit')
        ax1[1, 1].plot(xData, d1, 'r', label='doe1')
        ax1[1, 1].plot(xData, d2, 'c', label='doe2')
        ax1[1, 1].grid(True)
        ax1[1, 1].legend()
        ax1[1, 1].set_title(str([round(x, 2) for x in [A1, onset1, taur1, taud1]])
                            + '\n'
                            + str([round(x, 2) for x in [A2, onset2, taur2, taud2]]))
        ax1[1, 1].set_ylabel('Voltage (mV)')
        ax1[1, 1].set_xlabel('time (' + signal2Fit.t_start.units.dimensionality.string + ')')

        N = xData.shape[0]

        # https://www.researchgate.net/post/What_is_the_AIC_formula
        aic = N * np.log((fitError ** 2) / N) + 2 * 9

        bic = N * np.log((fitError ** 2) / N) + np.log(N) * 9

    aics_doe.append(aic)
    bics_doe.append(bic)

    fig1.tight_layout()
    fig1.canvas.draw()
    fig1.savefig(os.path.join(resDir, expName + '.eps'))

    raw_input()


fig2, ax2 = plt.subplots(figsize=(10, 8))
ax2.plot(range(len(expNames)), aics_alpha, 'bo', label='aics_alpha')
ax2.plot(range(len(expNames)), bics_alpha, 'ro', label='bics_alpha')
ax2.plot(range(len(expNames)), aics_doe, 'b^', label='aics_doe')
ax2.plot(range(len(expNames)), bics_doe, 'r^', label='bics_doe')

ax2.set_xticks(range(len(expNames)))
ax2.set_xticklabels(expNames, rotation=90)
ax2.legend()

fig2.tight_layout()
fig2.canvas.draw()
fig2.savefig(os.path.join(resDir, 'ModelComparision.eps'))