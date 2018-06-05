from GJEMS.ephys.rawDataAnalyse import RawDataAnalyser
from GJEMS.ephys.NEOFuncs import downSampleAnalogSignal, sliceAnalogSignal
from GJEMS.ephys.expIDLaborStateMap import expIDLaborStateMap
from GJEMS.ephys.doubleExpFitting import doubleExpFun, twoDoubleExps, accept_test1
from GJEMS.ephys.pdColumnNameMapCont import mdFN, fFN, paramIndex
from scipy.signal import lfilter, kaiserord, firwin
import numpy as np
from neo import AnalogSignal
from matplotlib import pyplot as plt
# plt.ion()
import os
import quantities as qu
import json
import operator
import pandas as pd


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
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting2Exp/')
toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')
scalesFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/scales.json')
tempParFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting2Exp/tempPars265All.json')

with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

with open(scalesFile, 'r') as fle:
    scales = json.load(fle)

traceLength = 1 * qu.s
# nPts = int(traceLength / Ts)

traceStart = 0 * qu.s
cutOff = 50
transitionWidth = 40
newSamplingRate = 400


types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [1 * qu.s, 1 * qu.s, 1 * qu.s]

fig1, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(14, 11.2))
fig, ax = plt.subplots(figsize=(14, 11.2))
fig2, ax2 = plt.subplots(figsize=(14, 11.2))

allData = pd.DataFrame(None, columns=fFN.values() + mdFN.values())
ax2_currentMax = 0 * qu.mV

for expName in expNames:

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

    resps = rda.getContResps(types=types)

    for freq, freqResps in resps.iteritems():
        print('Doing ' + str(freq) + 'Hz; All freqs ' + str(resps.keys()) + ' Hz')

        nrnFreqResDir = os.path.join(nrnResDir, expName + '-' + str(freq))

        normedFilteredSigs = []
        normedSigs = []
        t_starts = []
        ax2_currentMaxs = []

        for respInd, resp in enumerate(freqResps):

            print('{}/{}; trial starting at {}'.format(respInd + 1, len(freqResps), resp['DuringStimulus'].t_start))

            respNormedSigs = []
            respNormedSigsFiltered = []

            for typeInd, tpye in enumerate(types):

                temp = resp[tpye]
                if shouldBeIgnored(temp):
                    print('Trial{} {} ignored'.format(respInd + 1, tpye))
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

                delay, respFBFiltered = LPFilterKaiser(respFB, cutoff=cutOff, transitionWidth=transitionWidth, rippleDB=60)

                normedFilteredSigs.append(respFBFiltered)
                normedSigs.append(respFB)

                t_starts.append(resp['DuringStimulus'].t_start)


                downSamplingFactor = int(np.floor(respFBFiltered.sampling_rate / newSamplingRate))
                signal2FitFull = downSampleAnalogSignal(respFBFiltered, downSamplingFactor)
                signal2Fit = sliceAnalogSignal(signal2FitFull, traceStart, traceStart+traceLength)

                xDataT = signal2Fit.times
                xDataT.units = qu.ms
                xData = xDataT.magnitude
                yData = signal2Fit.magnitude

                outFile = os.path.join(nrnFreqResDir, '{}_Freq{}_Trial{}.json'.format(expName, freq, respInd))
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
                    else:
                        x = [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset]

                    xs[xInd] = x

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
                    xBest = [np.nan] * len(paramIndex)
                else:
                    xBest = xs_accepted[int(np.argsort(fs_accepted)[0])]

                [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest

                temp = {
                                'expID':        expName,
                                'freq':         freq,
                                'laborState':   expIDLaborStateMap(expName),
                                'trialName':    'Trial{:d}'.format(respInd),
                                'trialStart':   resp['DuringStimulus'].t_start.magnitude
                            }
                tempDFDict = {mdFN[k]: v for k, v in temp.iteritems()}
                assert len(tempDFDict) == len(mdFN), 'not all metadata covered in the final data file'

                for k, v in paramIndex.iteritems():
                    tempDFDict[fFN[k]] = xBest[v]

                allData = allData.append(pd.Series(tempDFDict), ignore_index=True)

                fitSignal = twoDoubleExps(xData, *xBest)
                d1 = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)
                d2 = -doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)
                fitError = np.linalg.norm(fitSignal - yData)

                print(d1.max(), -d2.min(), -d1.max() / d2.min())


                print(round(offset, 4))
                print(np.round([A1, onset1, taur1, taud1], 4).T)
                print(np.round([A2, onset2, taur2, taud2], 4).T)
                print(fitError)

                ax.clear()
                normedSig = respFB

                normedSigT = normedSig.times
                normedSigT.units = qu.ms
                normedFilteredSigT = respFBFiltered.times
                normedFilteredSigT.units = qu.ms
                ax.plot(normedSigT, normedSig, 'b', label='Sig')
                ax.plot(normedFilteredSigT, respFBFiltered, 'r', label='SigLP')
                ax.plot(xData, fitSignal, 'g', label='fit')
                ax.set_ylabel('Voltage (mV)')
                ax.set_xlabel('time (ms)')
                ax.legend()
                fig.canvas.draw()

                ax2_currentMax -= 1.25 * (normedSig.max() - normedSig.min())

                ax2.plot(normedSigT, normedSig + ax2_currentMax, 'b', label='Sig')
                ax2.plot(normedFilteredSigT, respFBFiltered + ax2_currentMax, 'r', label='SigLP')
                ax2.plot(xData, fitSignal + ax2_currentMax.magnitude, 'g', label='fit')

                # ax1.plot((centeredSignal.times - centeredSignal.t_start) * 1e3, centeredSignal, 'r', label='raw')
                ax1[0].plot(xData, yData, 'g', label='filtered')
                ax1[0].plot(xData, fitSignal, 'b', label='fit')
                ax1[0].text(xData[int(len(xData) * 0.15)], 0.75 * max(yData), 'SSE=' + str(round(fitError, 2)), fontsize=14)

                ax1[0].legend()
                ax1[0].grid(True)

                ax1[1].plot(xData, fitSignal, 'b', label='fit')
                ax1[1].plot(xData, d1, 'r', label='doe1')
                ax1[1].plot(xData, d2, 'c', label='doe2')
                ax1[1].grid(True)
                ax1[1].legend()
                ax1[1].set_ylabel('Voltage (mV)')
                ax1[1].set_xlabel('time (ms)')
                ax1[1].set_title(('doe1:', [round(x, 2) for x in [A1, onset1, taur1, taud1]], 'doe2:',
                                    [round(x, 2) for x in [A2, onset2, taur2, taud2]]))

                fig1.canvas.draw()
                fig1.savefig(os.path.join(nrnFreqResDir,
                                          expName + '_Freq' + str(freq) +'_Trial' + str(respInd) + '.eps'), dpi=300)
                fig1.savefig(os.path.join(nrnFreqResDir,
                                          expName + '_Freq' + str(freq) +'_Trial' + str(respInd) + '.jpg'), dpi=300)

                if enable_debugging:
                    ipdb.set_trace()
                # else:
                #     raw_input()


        ax2.set_ylabel('Voltage (mV)')
        ax2.set_xlabel('time (ms)')
        ax2.set_yticklabels([])
        fig2.canvas.draw()

        fig2.savefig(os.path.join(nrnResDir, expName + '_Freq' + str(freq) + 'allTrials.eps'), dpi=300)
        fig2.savefig(os.path.join(nrnResDir, expName + '_Freq' + str(freq) + 'allTrials.jpg'), dpi=300)

        if enable_debugging:
            ipdb.set_trace()
        # else:
        #     raw_input()

        ax2.clear()

allDataPivoted = allData.set_index(keys=[mdFN['expID'], mdFN['freq'], mdFN['trialName']])
allDataPivoted.to_excel(os.path.join(resDir, 'contStimAllData.xlsx'))









