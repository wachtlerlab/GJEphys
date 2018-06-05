from GJEphys.rawDataAnalyse import RawDataAnalyser
from GJEphys.pdColumnNameMapCont import mdFN, fFN, getXBestFromCurrentData, newFFN, newFFNUnits, extraParamFuncs
import pandas as pd
import os
import json
import quantities as qu
import operator
from neo import AnalogSignal, SpikeTrain
import numpy as np
from GJEMS.folderDefs import homeFolder




dirpath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting2Exp/')
toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')
allDataFile = os.path.join(resDir, 'contStimAllData.xlsx')

with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

traceLength = 1 * qu.s
# nPts = int(traceLength / Ts)

traceStart = 0 * qu.s

types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [3 * qu.s, 1 * qu.s, 1 * qu.s]





for k, epf in extraParamFuncs.iteritems():

    exec('from GJEMS.GJEphys.additionalEphysFuncs import {}'.format(epf))

allData = pd.read_excel(allDataFile, index_col=(0, 1, 2))
newData = pd.DataFrame(None, columns=newFFN.values(), index=allData.index.copy())

allRespData = {}
allSpikeData = {}

for expName, expDF in allData.groupby(level=0):
    expName = str(expName)
    print('Doing ' + expName)

    rda = RawDataAnalyser(expName, dirpath)

    resps = rda.getContResps(types=types)
    spikes = rda.getContSpikes(types=types)


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

    for freq, freqDF in expDF.groupby(level=1):
        print('Doing ' + str(freq) + 'Hz; All freqs ' + str(resps.keys()) + ' Hz')

        freqResps = resps[freq]
        normedSigs = []
        t_starts = []
        normedSigSpikes = []
        respsDone = 0

        for trial, trialDF in freqDF.groupby(level=2):
            trialInd = int(trial[5:])
            resp = freqResps[trialInd]
            print('{}/{}; trial starting at {}'.format(trialInd + 1, len(freqResps), resp['DuringStimulus'].t_start))

            respNormedSigs = []
            respNormedSigsFiltered = []
            respSpikes = []

            for typeInd, tpye in enumerate(types):

                temp = resp[tpye]
                if shouldBeIgnored(temp):
                    print('Trial{} {} ignored'.format(trialInd + 1, tpye))
                    break
                centeredSig = temp.magnitude
                typeLen = int(typeDurs[typeInd] * temp.sampling_rate)
                presSigScaled = np.zeros((typeLen,))
                sigLen = min(typeLen, temp.shape[0])


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
                respSpikes.append(spikes[freq][trialInd][tpye])

            if len(respNormedSigs) == 3:
                respsDone += 1
                respFBSignal = np.concatenate(respNormedSigs)
                baseLineSig = np.concatenate((respNormedSigs[0], respNormedSigs[2]))
                baseLine = np.median(baseLineSig)
                respFBSignal = respFBSignal - baseLine

                respFB = AnalogSignal(signal=respFBSignal,
                                      units=qu.mV,
                                      sampling_rate=temp.sampling_rate,
                                      t_start=-typeDurs[0],
                                      )

                normedSigs.append(respFB)


                allSpikeTimes = np.concatenate([sp.times.magnitude for sp in respSpikes]) * respSpikes[0].units
                allSpikeTimes -= resp['DuringStimulus'].t_start
                allSpikeTimes = allSpikeTimes[
                                np.logical_and(respFB.t_start <= allSpikeTimes, allSpikeTimes <= respFB.t_stop)]

                if len(allSpikeTimes):
                    spikeTimesUnits = allSpikeTimes[0].units
                else:
                    spikeTimesUnits = qu.s

                respSpikeTrain = SpikeTrain(times=allSpikeTimes,
                                            units=spikeTimesUnits,
                                            t_start=-typeDurs[0], t_stop=typeDurs[1] + typeDurs[2])

                currentData = trialDF.loc[expName, freq, trial]

                samplingRate = respFB.sampling_period
                samplingRate.units = qu.ms

                traceStart.units = qu.ms
                traceLength.units = qu.ms

                traceStartInMS = traceStart.magnitude
                traceLengthInMS = traceLength.magnitude

                xData = np.arange(traceStartInMS, traceStartInMS + traceLengthInMS, samplingRate.magnitude)

                basicArgs = [respFB, respSpikeTrain, getXBestFromCurrentData(currentData), xData]


                for k, v in newFFN.iteritems():

                    exec('temp = {}(*(basicArgs + [newFFNUnits[k]]))'.format(extraParamFuncs[k]))
                    try:
                        temp = float(temp)
                    except TypeError as e:
                        raise ValueError('Function {} does not return a float-convertable value'.format(
                                                                                                  extraParamFuncs[k]))
                    newData.loc[expName, freq, trial][v] = temp


allData = pd.concat([allData, newData], axis=1)
pre, post = allDataFile.split('.', 1)
allData.to_excel('{}_expanded.{}'.format(pre, post), sheet_name='Sheet1')



