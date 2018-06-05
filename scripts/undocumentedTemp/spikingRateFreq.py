from GJEMS.ephys.rawDataAnalyse import RawDataAnalyser
import numpy as np
from matplotlib import pyplot as plt
import quantities as qu
import ipdb

def getSpRate(sp):

    stDur = (sp.t_stop - sp.t_start).simplified.magnitude
    return sp.shape[0] / stDur

def getSpikesIn(sp, t_stop, t_start=0 * qu.ms):

    spikesIn = sp[(sp >= t_start + sp.t_start) & (sp <= sp.t_start + t_stop)]
    spikesIn.t_start = max(sp.t_start, t_start + sp.t_start)
    spikesIn.t_stop = min(sp.t_stop, t_stop + sp.t_start)

    return spikesIn



expNames = [
            # '130313-4Rh',
            # '130322-1LY',
            # '130326-2Rh',
            # '130408-1LY',
            # '130425-1Al',
            # '130501-2Rh',
            '130523-3LY',
            '130605-1LY',
            '130605-2LY',
            # '130705-1LY',
            # '140424-1LY',
            # '140701-1Al',
            '140813-3Al',
            '140930-1Al',
            ]

dirpath = 'NIXFiles'

fig1 = plt.figure()
plt.show(block=False)

cols = ['r', 'r', 'r', 'g', 'g', 'g', 'b', 'b', 'b']
ms = ['o', 's', '^', 'o', 's', '^', 'o', 's', '^']


ax1 = plt.subplot(221)
plt.title('Spike rate before the stimulus')
plt.xlabel('Frequency(Hz)')
plt.ylabel('Spiking Rate(spikes/s)')
plt.xlim([0, 410])
plt.ylim([0, 200])
plt.draw()

ax2 = plt.subplot(222)

plt.title('Spike rate during the stimulus [0,100]ms')
plt.xlabel('Frequency(Hz)')
plt.ylabel('Spiking Rate(spikes/s)')
plt.xlim([0, 410])
plt.ylim([0, 200])
plt.draw()

ax3 = plt.subplot(223)

plt.title('Spike rate during the stimulus [100, 1000]ms')
plt.xlabel('Frequency(Hz)')
plt.ylabel('Spiking Rate(spikes/s)')
plt.xlim([0, 410])
plt.ylim([0, 200])
plt.draw()


ax4 = plt.subplot(224)

plt.title('Spike rate after the stimulus')
plt.xlabel('Frequency(Hz)')
plt.ylabel('Spiking Rate(spikes/s)')
plt.xlim([0, 410])
plt.ylim([0, 200])
plt.draw()



allSpikesDuring = {}
allSpikesBefore = {}
allSpikesAfter = {}
for expInd, expName in enumerate(expNames):

    print('Reading ' + expName)

    rda = RawDataAnalyser(expName, dirpath)

    spikes = rda.getContSpikes()

    plt.subplot(221)
    plt.cla()
    plt.subplot(222)
    plt.cla()
    plt.subplot(223)
    plt.cla()
    plt.subplot(224)
    plt.cla()

    for freq, sp in spikes.iteritems():



        if freq in allSpikesDuring:

            allSpikesDuring[freq].extend(sp['DuringStimulus'])
            allSpikesBefore[freq].extend(sp['BeforeStimulus'])
            allSpikesAfter[freq].extend(sp['AfterStimulus'])


        else:

            allSpikesDuring[freq] = sp['DuringStimulus']
            allSpikesAfter[freq] = sp['AfterStimulus']
            allSpikesBefore[freq] = sp['BeforeStimulus']

        plt.subplot(221)
        plt.title('Spike rate before the stimulus')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spiking Rate(spikes/s)')
        plt.xlim([0, 410])
        plt.ylim([0, 200])

        for spt in sp['BeforeStimulus']:
            plt.plot(freq, getSpRate(spt), color=cols[expInd], marker=ms[expInd], ls='None', ms=5)


        for spt in sp['DuringStimulus']:

            plt.subplot(222)

            plt.title('Spike rate during the stimulus [0,50]ms')
            plt.xlabel('Frequency(Hz)')
            plt.ylabel('Spiking Rate(spikes/s)')
            plt.xlim([0, 410])
            plt.ylim([0, 200])

            spe = getSpikesIn(spt, t_stop=50 * qu.ms)
            plt.plot(freq, getSpRate(spe) + 5 * np.random.rand(1), color=cols[expInd], marker=ms[expInd], ls='None', ms=5)
            plt.subplot(223)

            plt.title('Spike rate during the stimulus [50, 1000]ms')
            plt.xlabel('Frequency(Hz)')
            plt.ylabel('Spiking Rate(spikes/s)')
            plt.xlim([0, 410])
            plt.ylim([0, 200])

            spi = getSpikesIn(spt, t_stop=1000 * qu.ms, t_start=50 * qu.ms)
            plt.plot(freq, getSpRate(spi), color=cols[expInd], marker=ms[expInd], ls='None', ms=5)


        plt.subplot(224)

        plt.title('Spike rate after the stimulus')
        plt.xlabel('Frequency(Hz)')
        plt.ylabel('Spiking Rate(spikes/s)')
        plt.xlim([0, 410])
        plt.ylim([0, 200])


        for spt in sp['AfterStimulus']:
            plt.plot(freq, getSpRate(spt), color=cols[expInd], marker=ms[expInd], ls='None', ms=5)

    plt.draw()
    ipdb.set_trace()


spikeRatesE = []
spikeRatesI = []
spikeRatesBefore = []
spikeRatesAfter = []






# for freq in allSpikesDuring.iterkeys():
#
#     # fig = plt.figure()
#     # plt.show(block=False)
#     #
#     # plt.suptitle(str(freq))
#
#     spRateE = []
#     spRateI = []
#
#     for sp in allSpikesDuring[freq]:
#         spE = getSpikesIn(sp, 100 * qu.ms)
#         spRateE.append(getSpRate(spE))
#         spI = getSpikesIn(sp, 1 * qu.s, 100 * qu.ms)
#         spRateI.append(getSpRate(spI))
#
#     spikeRatesE.append(spRateE)
#     spikeRatesI.append(spRateI)
#
#     spRateB = [getSpRate(sp) for sp in allSpikesBefore[freq]]
#     spikeRatesBefore.append(spRateB)
#
#     spRateA = [getSpRate(sp) for sp in allSpikesAfter[freq]]
#     spikeRatesAfter.append(spRateA)




# frMeanA = [np.mean(sp) for sp in spikeRatesAfter]
# frStdA = [np.std(sp) for sp in spikeRatesAfter]
# plt.errorbar(allSpikesDuring.keys(), frMeanA, yerr=frStdA, color='b', ls='None', marker='o')













