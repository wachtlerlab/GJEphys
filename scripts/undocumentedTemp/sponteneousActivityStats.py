import pandas as pd
import os
from GJEMS.folderDefs import homeFolder
from GJEphys.additionalEphysFuncs import spontAct3Sec
from GJEphys.pdColumnNameMapCont import mdFN, newFFN
from GJEphys.rawDataAnalyse import RawDataAnalyser
from GJEphys.KKHAXLParsing import getExpIDsByCategory
from GJEphys.expIDLaborStateMap import expIDLaborStateMap
import json
import operator
import quantities as qu
import numpy as np
from neo import SpikeTrain
import sys
from matplotlib import pyplot as plt
from GJEMS.viz.matplotlibRCParams import mplPars
import seaborn as sns

mplPars['axes.labelsize'] = "small"
sns.set(rc=mplPars, style="darkgrid")

types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
typeDurs = [3 * qu.s, 1 * qu.s, 1 * qu.s]

localAdditionalMDFN = {
    "prevFreq": "Frequency of preceeding stimulus (Hz)",
    "IStimI": "Time since \nprevious stimulus \n(s)"
}

def saveSpontActivityRates(outFile, toIgnoreFile, expNames, nixPath):

    with open(toIgnoreFile, 'r') as fle:
        toIgnore = json.load(fle)

    allData = pd.DataFrame()

    for expName in expNames:

        print('Doing ' + expName)

        if expName in toIgnore:
            intervalToIgnore = toIgnore[expName]
        else:
            intervalToIgnore = None

        def shouldBeIgnored(resp):

            if intervalToIgnore is None:
                return False
            else:
                respInIgnoreIntervals = [(x * qu.s <= resp.t_start <= y * qu.s) | (x * qu.s <= resp.t_stop <= y * qu.s)
                                            for x, y in intervalToIgnore]
                return reduce(operator.or_, respInIgnoreIntervals)

        rda = RawDataAnalyser(expName, nixPath)

        spikes = rda.getContSpikes(types=types)

        for freq, freqResps in spikes.iteritems():
            print('Doing ' + str(freq) + 'Hz; All freqs ' + str(spikes.keys()) + ' Hz')

            t_starts = []


            for stInd, sts in enumerate(freqResps):

                t_start = sts['DuringStimulus'].t_start
                print('{}/{}; trial starting at {}'.format(stInd + 1, len(freqResps),
                                                           t_start))

                respSpikes = []

                for typeInd, tpye in enumerate(types):

                    temp = sts[tpye]
                    if shouldBeIgnored(temp):
                        print('Trial{} {} ignored'.format(stInd + 1, tpye))
                        break

                    respSpikes.append(sts[tpye])

                if len(respSpikes) == 3:

                    t_starts.append(t_start)

                    allSpikeTimes = np.concatenate([sp.times.magnitude for sp in respSpikes]) * respSpikes[0].units
                    allSpikeTimes -= t_start
                    allSpikeTimes = allSpikeTimes[
                        np.logical_and(-typeDurs[0] <= allSpikeTimes, allSpikeTimes <= typeDurs[1] + typeDurs[2])]

                    if len(allSpikeTimes):
                        spikeTimesUnits = allSpikeTimes[0].units
                    else:
                        spikeTimesUnits = qu.s

                    respSpikeTrain = SpikeTrain(times=allSpikeTimes,
                                                units=spikeTimesUnits,
                                                t_start=-typeDurs[0], t_stop=typeDurs[1] + typeDurs[2])

                    temp = {
                        'expID': expName,
                        'freq': freq,
                        'laborState': expIDLaborStateMap(expName),
                        'trialName': 'Trial{:d}'.format(stInd),
                        'trialStart': t_start.magnitude
                    }
                    tempDFDict = {mdFN[k]: v for k, v in temp.iteritems()}
                    assert len(tempDFDict) == len(mdFN), 'not all metadata covered in the final data file'

                    tempDFDict[newFFN["spontFR3"]] = spontAct3Sec(None, respSpikeTrain, None, None)
                    allData = allData.append(pd.Series(tempDFDict), ignore_index=True)

    allDataPivoted = allData.set_index(keys=[mdFN['expID'], mdFN['freq'], mdFN['trialName']])
    allDataPivoted.to_excel(outFile)


def makeStatsImages(statsXL, outPDFFile):

    allDataOrig = pd.read_excel(statsXL, index_col=(0, 1, 2))
    allDataOrig.reset_index(inplace=True)

    figDir = outPDFFile[:-4]

    if not os.path.isdir(figDir):
        os.mkdir(figDir)

    allDataExpIDTimePivoted = allDataOrig.set_index(keys=[mdFN['expID'], mdFN['trialStart']])
    allDataExpIDTimePivoted.sort_index(level=(0, 1), inplace=True)

    allData = allDataExpIDTimePivoted.reset_index()

    allData = allData.set_index(keys=[mdFN['expID']])


    for expName, expDF in allData.groupby(level=0):

        temp = expDF.iloc[:-1, :][mdFN["freq"]].values
        temp1 = np.concatenate(([np.nan], temp))
        allData.loc[expName, localAdditionalMDFN["prevFreq"]] = temp1

        temp2 = np.diff(expDF[mdFN["trialStart"]].values)
        istimIs = np.concatenate(([np.nan], temp2))
        allData.loc[expName, localAdditionalMDFN["IStimI"]] = istimIs

    allData.reset_index(inplace=True)
    allData.to_excel("{}_extended.xlsx".format(statsXL[:-5]))

    for expName, expDF in allData.groupby(by=mdFN["expID"]):

        expFigFle = os.path.join(figDir, "{}.png".format(expName))
        fig, axs = plt.subplots(nrows=2, figsize=(12, 7.5))

        sns.violinplot(ax=axs[0], x=localAdditionalMDFN["prevFreq"],
                       y=localAdditionalMDFN["IStimI"], data=expDF)
        axs[0].set_ylim([-10, 50])
        sns.violinplot(ax=axs[1], x=localAdditionalMDFN["prevFreq"],
                       y=newFFN["spontFR3"], data=expDF)
        fig.suptitle(expName)

        fig.tight_layout()
        fig.savefig(expFigFle, dpi=150)











if __name__ == "__main__":

    assert len(sys.argv) == 2, "Invalid Usage. Please use as follows:\n " \
                               "python spontaneousActivityStats.py <option>\n" \
                               "Valid options are: saveXL, makeStatsImages"


    spike2Path = os.path.join(homeFolder, 'DataAndResults/GJEphys/spike2Files/')
    nixPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles')
    toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')
    excel = os.path.join(homeFolder, 'DataAndResults/GJEphys/spike2Files/neuron_database_20150720ha_ver6_modAK.xlsx')
    excelSheet = 'Kai-san final report150803'
    resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Results/')

    filterExpIDs = os.path.join(homeFolder, 'DataAndResults/GJEphys/expIDFilters/Ai2017_DL-Int-1.xlsx')

    filterExpIDsDF = pd.read_excel(filterExpIDs)

    expIDsByCat = getExpIDsByCategory(excel, excelSheet, ["DL-Int-1"], spike2Path)

    expNames = filterExpIDsDF.loc[:, mdFN["expID"]]

    filterStub = os.path.split(filterExpIDs)[1][:-5]
    filterOutDir = os.path.join(resDir, filterStub)
    if not os.path.isdir(filterOutDir):
        os.mkdir(filterOutDir)
    outFile = os.path.join(filterOutDir, "spontAct3Sec.xlsx")
    outPDFFile = os.path.join(filterOutDir, "spontAct3SecVsPrevFreqs.pdf")

    if sys.argv[1] == "saveXL":

        saveSpontActivityRates(outFile=outFile,
                               toIgnoreFile=toIgnoreFile,
                               nixPath=nixPath,
                               expNames=expNames)

    elif sys.argv[1] == "makeStatsImages":

        makeStatsImages(outFile, outPDFFile)

    else:

        raise(ValueError("Invalid input {} for <option>".format(sys.argv[1])))