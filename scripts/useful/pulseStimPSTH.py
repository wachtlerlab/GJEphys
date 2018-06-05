'''
Description: This script is used to generate raster and PTSH plots of pulse train stimuli. For every <Experiment ID>
specified, subplots are made, one per valid combination of (Pulse Duration, Pulse Interval), containing the raster and
PSTH plots of responses of all trials. The resulting figure is saved into PNG and EPS files.

The inputs can be specified by modifying the following variables in the code below:

1. expNames: list of <Experiment ID>s to use.
2. dirPath: string containing the path the directory containing processed NIX Files
3. outputDir: string containing a valid path in the file system where the output PNGs and EPSs, one per <Experiment ID>
specified, is created.

Usage:
python <path to parent directory>/pulseStimPSTH.py


'''

import os
import nixio as nix
homeFolder = os.path.expanduser('~')
from GJEphys.neoNIXIO import property2qu
import quantities as qu
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
# plt.ion()
import json

def getSecPulseProps(section, units):

    if 'PulseDuration' in section:
        pulseDur = property2qu(section.props['PulseDuration'])
        pulseDur.units = units
    else:
        raise('Pulse Duration not found in ' + section.name)

    if 'PulseInterval' in section:
        pulseInt = property2qu(section.props['PulseInterval'])
        pulseInt.units = units
    else:
        raise('Pulse Interval not found in ' + section.name)

    return float(pulseDur.magnitude), float(pulseInt.magnitude)

def getTagPositionQu(tag, units):

    pos = qu.Quantity(tag.position, units=tag.units[0])
    pos.units = units
    return pos

def getMultitagPositions(multitag, units):

    pos = qu.Quantity(multitag.positions, units=multitag.units[0])
    pos.units = units

    return pos



expNames = [
            # '130313-4Rh',
            # '130322-1LY',
            # '130326-2Rh',
            # '130408-1LY',
            # '130425-1Al',
            '130501-2Rh',
            # '130523-3LY',
            # '130605-1LY',
            # '130605-2LY',
            '130705-1LY',
            # '140424-1LY',
            # '140701-1Al',
            # '140813-3Al',
            # '140930-1Al',
            # '140917-1Al',
            # '141030-1Al',
            #
            '130819-2LY',
            # '130822-3Al'
            # '130415-2Rh'

            # '130318-3LY',
            # '130425-3Al',
            # '130517-1Al',
            # '130704-1LY',
]

dirPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles')
outputDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Results/PulseStimUniModality102030b50')

if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

# pulseParams = [
#                 (10.0, 50.0),
#                 (16.0, 100.0),
#                 (16.0, 50.0),
#                 (16.0, 33.0),
#                 (20.0, 50.0),
#                 (27.0, 50.0),
#                 (30.0, 50.0),
#                 (30.0, 150.0),
#                 (30.0, 100.0),
#                 (30.0, 75.0),
#                 (53.0, 100.0)
#                 ]

pulseParams = [(10.0, 50.0), (20.0, 50.0), (30.0, 50.0)]


spikeLatMeans = {k:[] for k in pulseParams}
spikeLatStds = {k:[] for k in pulseParams}
spikeLatExpInds = {k:[] for k in pulseParams}

durations = sorted(list(set([x[0] for x in pulseParams])))
intervals = sorted(list(set([x[1] for x in pulseParams])))

pulseParUnits = qu.ms
spikeLatencyUnits = qu.ms

fig0, ax0 = plt.subplots(1, 3, figsize=(10, 7.5))


pulseExpInds = []
for expInd, expName in enumerate(expNames):

    spikeLatencies = {k:[] for k in pulseParams}
    trialIndMax = 0

    print('Doing ' + expName)
    fName = os.path.join(dirPath, expName + '.h5')

    nixFile = nix.File.open(fName, nix.FileMode.ReadOnly)

    if 'PulseStimulii' in nixFile.sections['VibrationStimulii-Processed'].sections:


        for tag in nixFile.blocks['RawDataTraces'].tags:

            if tag.type[:7] == 'OnPulse':



                stimStartQ = getTagPositionQu(tag, spikeLatencyUnits)

                pulsePars = getSecPulseProps(tag.metadata.parent, pulseParUnits)

                mtagName = 'MultiTag' + tag.name[3:]


                if pulsePars in pulseParams:

                    trialInd = int(tag.metadata.name[5:])
                    pulseInd = int(tag.type[7:])
                    trialIndMax = max(trialInd, trialIndMax)
                    ax = ax0[pulseParams.index(pulsePars)]

                    if expInd not in pulseExpInds:
                            pulseExpInds.append(expInd)
                            [ax.clear() for ax in ax0]

                    if mtagName in nixFile.blocks['RawDataTraces'].multi_tags:
                        spikeTimesQ = getMultitagPositions(
                                        nixFile.blocks['RawDataTraces'].multi_tags[mtagName],
                                        spikeLatencyUnits)



                        spikeLatenciesT = (spikeTimesQ - stimStartQ).magnitude.tolist()

                        spikeLatencies[pulsePars].extend(spikeLatenciesT)





                        ax.plot(spikeLatenciesT,
                                [3 * trialInd + 0.1 * pulseInd] * len(spikeLatenciesT), 'r^', ms=5)
                        ax.plot([0, max(durations)], [3 * trialInd + 0.1 * pulseInd] * 2, 'r-', lw=0.5)
                        ax.set_xlabel('Time in ' + spikeLatencyUnits.dimensionality.string)

                    else:

                        ax.plot([0, max(durations)], [3 * trialInd + 0.1 * pulseInd] * 2, 'g-', lw=0.5)

        expJSON = {}

        for k, v in spikeLatencies.iteritems():
            if v:
                spikeLatMeans[k].append(float(np.mean(v)))
                spikeLatStds[k].append(float(np.std(v)))
                spikeLatExpInds[k].append(expInd)

                ax = ax0[pulseParams.index(k)]
                pulseParamUnitStr = pulseParUnits.dimensionality.string
                ax.set_title({'dur': str(k[0]) + pulseParamUnitStr, 'int': str(k[1]) + pulseParamUnitStr})
                hist, bins = np.histogram(v, bins=np.arange(max(durations)))
                spikeYMax = 3 * trialIndMax + 1
                histYExtent = round(0.5 * spikeYMax)
                ax.bar(bins[:-1], hist * histYExtent /  max(hist), bottom=- histYExtent)
                ax.set_yticks(np.linspace(-histYExtent, 0, max(hist) + 1))
                ax.set_yticklabels([str(x) for x in range(max(hist) + 1)])



        [ax.set_xlim(0, max(durations)) for ax in ax0]
        # fig0.tight_layout()
        # fig0.canvas.draw()
        fig0.savefig(os.path.join(outputDir, expName + '_PulseStims.eps'))
        fig0.savefig(os.path.join(outputDir, expName + '_PulseStims.png'), dpi=200)





