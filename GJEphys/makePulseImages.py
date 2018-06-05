import sys
import os
from matplotlib import pyplot as plt
import quantities as qu
from GJEphys.rawDataAnalyse import RawDataAnalyser
import json
import seaborn as sns
import numpy as np

assert len(sys.argv) == 2, 'Improper Usage! Please use as follows:\npython {} <json param file>'.format(sys.argv[1])

with open(sys.argv[1]) as fle:

    pars = json.load(fle)
    NIXPath = pars['NIXPath']
    expName = pars['expName']
    catResDir = pars['catResDir']
    downSampleFactor = pars['downSampleFactor']
    type2Color = pars['type2Color']
    mplPars = pars['mplPars']

sns.set(rc=mplPars)

rda = RawDataAnalyser(expName, NIXPath)

pulseResps = rda.getPulseResps()
pulseSpikes = rda.getPulseSpikes()

if pulseResps and pulseSpikes:

    expDir = os.path.join(catResDir, expName)

    if not os.path.isdir(expDir):
        os.mkdir(expDir)

    allPPs = pulseResps.keys()
    print('Found these pulse Pars:{}'.format(allPPs))

    allPulseDurs = list(sorted(set([x[0] for x in allPPs])))
    allPulseInts = list(sorted(set([x[1] for x in allPPs])))
    temp = list(set([x[2] for x in allPPs]))
    assert len(temp) == 1, "more than one base pulse frequency found. Something's wrong!"
    basePulseFreq = temp[0]


    with sns.axes_style('darkgrid'):
        overallFig, overallAxs = plt.subplots(nrows=len(allPulseDurs), ncols=len(allPulseInts), figsize=(7, 5.6))

        overallAxs = np.array(overallAxs).reshape((len(allPulseDurs), len(allPulseInts)))

        overallFig.text(0.00, 0.5, 'Pulse Duration (ms)', va='center', rotation='vertical')
        overallFig.text(0.5, 0.98, 'Pulse Interval (ms)', va='center', ha='center')
        overallFig.text(0.97, 0.5, 'Trial No', va='center', rotation='vertical')
        overallFig.text(0.5, 0.02, 'Pulse Number in order of application', va='center', ha='center')

        for rowInd in range(overallAxs.shape[0]):
            for colInd in range(overallAxs.shape[1]):

                dur, inter = allPulseDurs[rowInd], allPulseInts[colInd]
                ax = overallAxs[rowInd][colInd]

                ax.set_ylim([-1, 11])
                ax.set_xlim(-inter, 11 * inter)
                ax.set_xticks(np.arange(0, 10 * inter, inter))
                ax.set_yticks(range(11))


                if colInd < overallAxs.shape[1] - 1:
                    ax.set_yticklabels([])
                else:
                    ax.yaxis.tick_right()
                    ax.set_yticklabels([str(x) if x % 2 else '' for x in range(1, 11)])

                if colInd == 0:

                    ax.set_ylabel(dur)

                if rowInd == 0:

                    ax.set_title(inter)

                if rowInd < overallAxs.shape[0] - 1:

                    ax.set_xticklabels([])
                else:
                    ax.set_xticklabels([str(x) if x % 2 else '' for x in range(1, 11)])

                if (dur, inter, basePulseFreq) not in allPPs:
                    ax.text(0.5, 0.5, 'Not \napplied', va='center', ha='center', color='r', transform=ax.transAxes)





    for pp, ppResps in pulseResps.iteritems():

        print('Doing {}'.format(pp))

        ppDir = os.path.join(expDir, str(pp))

        if not os.path.isdir(ppDir):
            os.mkdir(ppDir)

        ppAx = overallAxs[allPulseDurs.index(pp[0])][allPulseInts.index(pp[1])]

        for trialInd, trialResps in enumerate(ppResps):

            trialStart = trialResps['DuringStimulus'].t_start
            with sns.axes_style('darkgrid'):
                fig, axs = plt.subplots(nrows=2, figsize=(7, 5.6))
                fig.suptitle(r'Duration={}, Interval={}, Trial\#={}'.format(pp[0], pp[1], trialInd + 1))

                for typ, typResps in trialResps.iteritems():

                    if typ in ['OnPulses', 'OffPulses', 'BeforeStimulus', 'AfterStimulus']:

                        typSpikes = pulseSpikes[pp][trialInd][typ]

                        if type(typResps) is list:

                            for respInd, resp in enumerate(typResps):

                                respSpikes = typSpikes[respInd]

                                axs[0].plot(resp.times[::downSampleFactor], resp[::downSampleFactor],
                                            color=type2Color[typ], marker=None, ls='-')
                                axs[1].plot(resp.times[::downSampleFactor], resp[::downSampleFactor],
                                            color=type2Color[typ], marker=None, ls='-')
                                axs[0].plot(respSpikes.times, [resp.max()] * respSpikes.size, '^k')
                                axs[1].plot(respSpikes.times, [resp.max()] * respSpikes.size, '^k')

                                spTimes = qu.Quantity(respSpikes.times) - trialStart
                                spTimes.units = qu.ms
                                ppAx.plot(spTimes, [trialInd] * respSpikes.size,
                                          color=type2Color[typ], marker='^', ls='None')
                        else:

                            axs[0].plot(typResps.times[::downSampleFactor], typResps[::downSampleFactor],
                                        color=type2Color[typ], marker=None, ls='-')
                            axs[0].plot(typSpikes.times, [typResps.max()] * typSpikes.size, '^k')

                for k, v in type2Color.iteritems():

                    axs[0].plot([], [], color=v, label=k, marker=None, ls='-')

                axs[1].set_xlabel('time(s)')
                axs[1].set_ylabel('Membrane \nPotential (mV)')
                # axs[0].legend(loc=(1.15, 1), ncol=1)
                fig.tight_layout()
                fig.savefig(os.path.join(ppDir, 'Trial{}.png'.format(trialInd + 1)), dpi=300)

            plt.close(fig.number)
            del fig

    del rda
    del pulseResps
    del pulseSpikes

    overallFig.tight_layout()
    overallFig.savefig(os.path.join(catResDir, '{}Overview.png'.format(expName)), dpi=300)

    plt.close(overallFig.number)
    del overallFig

