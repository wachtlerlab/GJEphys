import nixio as nix
import os
from GJEMS.folderDefs import homeFolder
import numpy as np
import quantities as qu
from GJEphys.neoNIXIO import property2qu, multiTag2SpikeTrain, simpleFloat
from matplotlib import pyplot as plt
import seaborn as sns
from GJEphys.matplotlibRCParams import mplPars
from GJEphys.folderDefs import allExpIDsByCat
from GJEphys.pdColumnNameMapPulse import mdFN, fFN
import pandas as pd

dirPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles')
sns.set(rc=mplPars)


def sectionDurIntervalTuple(sec, units):
    if 'PulseDuration' in sec.props:
        duration = property2qu(sec.props['PulseDuration'])[0]
    else:
        return None

    if 'PulseInterval' in sec.props:
        interval = property2qu(sec.props['PulseInterval'])[0]
    else:
        return None

    duration.units = units
    interval.units = units

    return float(duration), float(interval)

def isPulseTagF(multiTag):
    tagType = multiTag.type
    if tagType.startswith("OnPulse"):
        respType = "OnPulse"
        pulseNumber = int(tagType[7:])
        return True, respType, pulseNumber
    elif tagType.startswith("OffPulse"):
        respType = "OffPulse"
        pulseNumber = int(tagType[8:])
        return True, respType, pulseNumber
    else:
        return False, None, None

def getMultiTagStartEnd(multiTag, block):

    tagName = 'Tag{}'.format(multiTag.name[8:])
    tag = block.tags[tagName]
    stimStart = qu.Quantity(tag.position, units=tag.units[0])
    stimDur = qu.Quantity(tag.extent, units=tag.units[0])
    stimEnd = stimStart + stimDur

    return stimStart, stimEnd


def getMTSpikeTimes(multiTag, block):

    stStart, stEnd = getMultiTagStartEnd(multiTag, block)
    return multiTag2SpikeTrain(multiTag, stStart, stEnd)



expIDsByCat = allExpIDsByCat()
# expIDsByCat = {"DL-Int-1": ["130313-4Rh"]}
# expIDsByCat = {"DL-Int-1": ["130408-1LY"]}

cats = ["DL-Int-1"]


stimParams = [
                     (16.0, 33.0),
                     (30.0, 50.0),
                     (10.0, 50.0),
                     (20.0, 50.0),
                     (27.0, 50.0),
                     (16.0, 50.0),
                     (30.0, 75.0),
                     (30.0, 100.0),
                     (53.0, 100.0),
                     (16.0, 100.0),
             ]

stimParamUnits = qu.ms
outDir = os.path.join(homeFolder, "DataAndResults", "ephys", "Results", "PulseParamTuning")
if not os.path.isdir(outDir):
    os.makedirs(outDir)


for cat in cats:
    print("Doing {}".format(cat))

    catDir = os.path.join(outDir, cat)

    if not os.path.isdir(catDir):
        os.mkdir(catDir)

    expNames = expIDsByCat[cat]
    stimRespsData = pd.DataFrame()

    for expInd, expName in enumerate(expNames):

        print('Doing ' + expName)
        fName = os.path.join(dirPath, expName + '.h5')
        nixFile = nix.File.open(fName, nix.FileMode.ReadOnly)

        if 'PulseStimulii' in nixFile.sections['VibrationStimulii-Processed'].sections:
            for multiTag in nixFile.blocks['RawDataTraces'].multi_tags:

                isPulseTag, respType, pulseNumber = isPulseTagF(multiTag)

                if isPulseTag:
                    stimKey = sectionDurIntervalTuple(multiTag.metadata.parent, stimParamUnits)

                    if stimKey in stimParams:

                        tempDict = {}
                        tempDict[mdFN["respType"]] = respType
                        tempDict[mdFN["expID"]] = expName
                        tempDict[mdFN["cat"]] = cat
                        stimKey = sectionDurIntervalTuple(multiTag.metadata.parent, stimParamUnits)
                        tempDict[mdFN["pulseDur"]] = stimKey[0]
                        tempDict[mdFN["pulseInt"]] = stimKey[1]
                        tempDict[mdFN["pulseNumber"]] = pulseNumber

                        spikeTrain = getMTSpikeTimes(multiTag, nixFile.blocks['RawDataTraces'])
                        relSpikeTimes = qu.Quantity(spikeTrain.times) - spikeTrain.t_start
                        relSpikeTimesMS = simpleFloat(relSpikeTimes / qu.ms)
                        tempDict[fFN["totalSpikeNumber"]] = spikeTrain.times.shape[0]
                        if len(relSpikeTimesMS):
                            tempDict[fFN["firstSpikeLat"]] = relSpikeTimesMS[0]
                        else:
                            tempDict[fFN["firstSpikeLat"]] = np.nan
                        tempDict[fFN["allSpikes"]] = relSpikeTimesMS

                        stimRespsData = stimRespsData.append(tempDict, ignore_index=True)


    onStimResps = stimRespsData[stimRespsData[mdFN["respType"]] == "OnPulse"]
    fig1, ax1 = plt.subplots(figsize=(14, 11.2))
    fig2, ax2 = plt.subplots(figsize=(14, 11.2))
    sns.factorplot(y=fFN["firstSpikeLat"], x=mdFN["pulseDur"], hue=mdFN["pulseInt"], data=onStimResps,
                   kind="violin", ax=ax1)
    sns.factorplot(y=fFN["totalSpikeNumber"], x=mdFN["pulseDur"], hue=mdFN["pulseInt"], data=onStimResps,
                   kind="violin", ax=ax2)
    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig(os.path.join(catDir, "firstSpikeLat.png"), dpi=150)
    fig2.savefig(os.path.join(catDir, "totalSpikeNumber.png"), dpi=150)
    

    stimRespsData.to_excel(os.path.join(catDir, "AllPulseRespData.xlsx"))
    onStimResps.to_excel(os.path.join(catDir, "OnPulseRespData.xlsx"))