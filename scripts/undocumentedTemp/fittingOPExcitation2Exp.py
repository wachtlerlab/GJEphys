from GJEMS.ephys.rawDataAnalyse import RawDataAnalyser
from GJEMS.ephys.NEOFuncs import downSampleAnalogSignal, sliceAnalogSignal
from GJEMS.ephys.doubleExpFitting import runManyDoubleExpFittings
from scipy.signal import lfilter, kaiserord, firwin
import numpy as np
from neo import AnalogSignal
import os
import quantities as qu
import json
import operator
import shutil


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

if not os.path.isdir(resDir):
    os.mkdir(resDir)


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


nIterations = 1000


types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [1 * qu.s, 1 * qu.s, 1 * qu.s]


for expName in expNames:

    print('Doing ' + expName)

    nrnResDir = os.path.join(resDir, expName)

    if os.path.exists(nrnResDir):
        shutil.rmtree(nrnResDir)

    os.mkdir(nrnResDir)

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

        print('Doing {}Hz; All freqs {}Hz'.format(freq, resps.keys()))

        nrnFreqResDir = os.path.join(nrnResDir, expName + '-' + str(freq))

        if os.path.exists(nrnFreqResDir):
            shutil.rmtree(nrnFreqResDir)

        os.mkdir(nrnFreqResDir)

        normedFilteredSigs = []
        normedSigs = []

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
                delay, respFBFiltered = LPFilterKaiser(respFB, cutoff=cutOff,
                                                       transitionWidth=transitionWidth, rippleDB=60)

                downSamplingFactor = int(np.floor(respFBFiltered.sampling_rate / newSamplingRate))
                signal2FitFull = downSampleAnalogSignal(respFBFiltered, downSamplingFactor)
                signal2Fit = sliceAnalogSignal(signal2FitFull, traceStart, traceStart + traceLength)

                xDataT = signal2Fit.times
                xDataT.units = qu.ms
                xData = xDataT.magnitude
                yData = signal2Fit.magnitude

                outFile = os.path.join(nrnFreqResDir,
                                       '{}_Freq{}_Trial{}.json'.format(expName, freq, respInd))

                runManyDoubleExpFittings(xData, yData, nIterations, outFile, tempParFile)





















