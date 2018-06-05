from GJEphys.rawDataAnalyse import RawDataAnalyser
from GJEphys.folderDefs import NIXPath
import numpy as np
from matplotlib import pyplot as plt


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
            # # # '140917-1Al',
            # # # '141030-1Al',
            ]

fnew = [1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0]


allResps = {}
freqs = set()

for expName in expNames:

    print('Doing ' + expName)

    rda = RawDataAnalyser(expName, NIXPath)

    resps = rda.getContResps(types=['DuringStimulus'])

    allResps[expName] = resps

    freqs = freqs | (resps.viewkeys() - freqs)


avgTracesAll = [{}, {}]
traceLength = 1
Ts = 4.8e-5
nPts = int(0.9 / Ts)

for expName, expResp in allResps.iteritems():

    avgTraces = avgTracesAll[fnew[expNames.index(expName)]]

    for freq, freqResps in expResp.iteritems():

        for trial in freqResps:

            analogSignal = trial['DuringStimulus']

            if analogSignal.shape[0] >= nPts:

                if freq in avgTraces:

                    avgTraces[freq] = np.vstack((avgTraces[freq], analogSignal.magnitude[:nPts]))

                else:

                    avgTraces[freq] = analogSignal.magnitude[:nPts]

t = np.arange(nPts) * Ts
plt.ion()
for freq in freqs:

    fig = plt.figure()

    plt.suptitle('Frequency= ' + str(freq))

    if freq == 265:
        plt.subplot(2, 1, 1)

        mn = avgTracesAll[0][freq].mean(axis=0)
        std = avgTracesAll[0][freq].std(axis=0)
        plt.plot(t, mn, 'r')

        for ind in range(avgTracesAll[0][freq].shape[0]):
            plt.plot(t, avgTracesAll[0][freq][ind, :], color=[1, 0, 0, 0.1], marker='None', ls='-')
        plt.xlabel('time(s)')
        plt.ylabel('Voltage(mV)')
        plt.title('Newly Emerged(n=' + str(avgTracesAll[0][freq].shape[0]) + ')')

    if freq == 265:
        plt.subplot(2, 1, 2)

    if freq in avgTracesAll[1]:
        print([freq, avgTracesAll[1][freq].shape])
        mn = avgTracesAll[1][freq].mean(axis=0)
        std = avgTracesAll[1][freq].std(axis=0)
        plt.plot(t, mn, 'b')
        for ind in range(avgTracesAll[1][freq].shape[0]):
            plt.plot(t, avgTracesAll[1][freq][ind, :], color=[0, 0, 1, 0.1], marker='None', ls='-')
        plt.xlabel('time(s)')
        plt.ylabel('Voltage(mV)')
        plt.title('Forager(n=' + str(avgTracesAll[1][freq].shape[0]) + ')')













