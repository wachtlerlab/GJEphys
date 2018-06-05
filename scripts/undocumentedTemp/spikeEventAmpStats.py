from GJEphys.rawDataProcess import RawDataProcessor
from GJEphys.spikeDetection import SpikeDetector, getLocalMaximaThresholded, getGMMLastClusterCut, getLocalMaxima
import matplotlib.pyplot as plt
import quantities as qu
import numpy as np
plt.ion()

def thresholdAndDetect(filteredSignal, minSpikeWidth, maxSpikeWidth, thres):

        maxSpikeWidthSamples = (maxSpikeWidth * qu.ms / filteredSignal.sampling_period).simplified.magnitude

        minSpikeWidthSamples = (minSpikeWidth * qu.ms / filteredSignal.sampling_period).simplified.magnitude

        thresholded = np.int8(filteredSignal > thres)

        thresholdedDiff = np.diff(thresholded)

        onsets = np.where(thresholdedDiff > 0)[0] + 1

        offsets = np.where(thresholdedDiff < 0)[0] + 1

        spikeInds = []
        spikeHeights = []

        for onset, offset in zip(onsets, offsets):

            durAboveThesh = offset - onset
            if (durAboveThesh < maxSpikeWidthSamples) and (durAboveThesh >= minSpikeWidthSamples):

                spike = filteredSignal[onset:offset + 1]
                spikeHeight = spike.max().magnitude
                spikeHeights.append(spikeHeight)
                spikeInd = (spike.argmax() + onset)
                spikeInds.append(spikeInd)


        spikeTimes = filteredSignal.t_start + np.array(spikeInds) * filteredSignal.sampling_period
        return spikeTimes, np.array(spikeHeights) * filteredSignal.units

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

# assert len(expName) == 1

for expName in expNames:
    fig, ax = plt.subplots()
    print('Processing ' + expName)
    see = RawDataProcessor(expName=expName, dirpath='scripts/NIXFiles', readOnly=True)
    sd = SpikeDetector(see.voltageSignal)
    sd.filterButterworth()

    spikeTimes, spikeHeights = thresholdAndDetect(sd.filteredSignal,
                                                  minSpikeWidth=0.5, maxSpikeWidth=5, thres= sd.getNoiseStd())

    # maxInds, maxs = getLocalMaximaThresholded(sd.analogSignal.magnitude, 2 * sd.getNoiseStd())
    # maxInds, maxs = getLocalMaxima(sd.filteredSignal.magnitude)

    clusterCut = getGMMLastClusterCut(spikeHeights.magnitude, nComponents=2)
    # spikeTimes, spikeHeights = thresholdAndDetect(sd.filteredSignal,
    #                                               minSpikeWidth=0.5, maxSpikeWidth=5, thres=3 * sd.getNoiseStd())
    # see.close()
    # hist, bins = np.histogram(spikeHeights, bins=20)
    #
    # ax.cla()
    # ax.bar(bins[:-1], hist, width=bins[1] - bins[0])
    # plt.title(sd.getNoiseStd())
    # plt.draw()
    raw_input()
    plt.close('all')
