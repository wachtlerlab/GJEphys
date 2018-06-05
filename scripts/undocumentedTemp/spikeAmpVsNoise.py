from GJEphys.rawDataAnalyse import RawDataAnalyser
import numpy as np
from matplotlib import pyplot as plt
plt.ion()
import os
import quantities as qu
import json
import operator

homeFolder = os.path.expanduser('~')

mplPars = { 'text.usetex'       :    True,
            'axes.labelsize'    :   'large',
            'font.family'       :   'serif',
            'font.sans-serif'   :   'computer modern roman',
            }

for a, b in mplPars.items():
            plt.rcParams[a] = b

# ----------------------------------------------------------------------------------------------------------------------

def getSpikeAmps(resp, spikeTimes):

    spikeInds = map(int, (spikeTimes - resp.t_start) * resp.sampling_rate)
    return resp[spikeInds]

def shouldBeIgnored(resp, intervalToIgnore):

        if intervalToIgnore is None:
            return False
        else:
            respInIgnoreIntervals = [(x * qu.s <= resp.t_start <= y * qu.s) | (x * qu.s <= resp.t_stop <= y * qu.s)
                                        for x, y in intervalToIgnore]
            return reduce(operator.or_, respInIgnoreIntervals)

def getNoiseVar(resp):

    temp = resp - np.median(resp)
    return np.median(np.abs(temp)) / 0.6745


def fitLine(x, y):

    assert len(x) == len(y)

    matA = np.ones([len(x), 2])
    matB = np.asarray(y).T

    matA[:, 0] = np.asarray(x).T

    m, c = np.linalg.lstsq(matA, matB)[0]

    return m, c

# ----------------------------------------------------------------------------------------------------------------------

expNames = [
            '130313-4Rh',
            '130322-1LY',
            '130326-2Rh',
            '130408-1LY',
            '130425-1Al',
            '130501-2Rh',
            '130705-1LY',
            '140424-1LY',

            '130523-3LY',
            '130605-1LY',
            '130605-2LY',
            '140813-3Al',
            '140930-1Al',
            '140917-1Al',
            '141030-1Al',

            # '140701-1Al',

            # '130318-3LY',
            # '130425-3Al',
            # '130517-1Al',
            # '130704-1LY',
            ]



dirpath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
saveDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Results/responses265')
toIgnoreFile = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/toIgnore.json')
if not os.path.isdir(saveDir):
    os.mkdir(saveDir)

fig1, ax1 = plt.subplots(figsize=(10, 8))

axs = [ax1]
figs = [fig1]

Ts = 4.8e-5
respInitLen = 0.2
nPtsRespInitLen = int(respInitLen / Ts)

befSigLen = 0.1
nPtsBefLen = int(befSigLen / Ts)


t = np.arange(-nPtsBefLen, nPtsRespInitLen) * Ts


with open(toIgnoreFile, 'r') as fle:
    toIgnore = json.load(fle)

types=['DuringStimulus']
freqs = [10, 50, 100, 200, 265, 300, 400]

snrs = []
spikeAmpsAll = []
noiseAmpsAll = []

for expName in expNames:

    for ax in axs:
        ax.clear()

    print('Doing ' + expName)

    intervalToIgnore = None
    if expName in toIgnore:
        intervalToIgnore = toIgnore[expName]

    rda = RawDataAnalyser(expName, dirpath)

    resps = rda.getContResps(types=types)
    spikes = rda.getContSpikes(types=types)

    spikeAmps = []
    noiseAmps = []

    for freqInd, freq in enumerate(freqs):
        if freq in resps:
            for typeInd, tpye in enumerate(types):

                for respInd, respAll in enumerate(resps[freq]):

                    resp = respAll[tpye]

                    if not shouldBeIgnored(resp, intervalToIgnore):

                        respSpikeAmps = getSpikeAmps(resp, spikes[freq][respInd][tpye].times) \
                                        - np.median(resp.magnitude) * resp.units
                        noiseAmp = getNoiseVar(resp.magnitude)

                        spikeAmps.extend(respSpikeAmps)
                        noiseAmps.append(noiseAmp)
    snr = np.mean(spikeAmps) / np.mean(noiseAmps)
    snrs.append(snr)
    spikeAmpsAll.append(spikeAmps)
    noiseAmpsAll.append(noiseAmps)

ax1.plot(snrs, 'bo')
ax1.set_xticks(range(len(snrs)))
ax1.set_xticklabels(expNames, rotation=90)
ax1.set_ylabel('SNR')

fig1.tight_layout()
fig1.canvas.draw()




