'''
Description: This script is used to make and save individual figures of responses to continuous stimuli
as well as overview figures. The input parameters required for generating such plots are provided in a json file which is
supplied as a commandline argument. The json file must contain a dictionary with the following entries:
"NIXPath": <string containing the file system path of a directory containing processed NIX files>
"expName": <string containing a valid Experiment ID>
"freqs": <list of floats, indicating the frequencies of stimulus to use, if they are present in the experiment,
"catResDir": <string, indicating the file system path of an existing folder, within which generated images will be
hierarchically organized. A folder with <Experiment ID> will be created within it, and with the <Experiment ID> folder,
folders will be created, one for each present frequency>
"downSampleFactor": float, factor by which traces are temporally downsampled before plotting.
"type2Color": dict, must have the keys "BeforeStimulus", "DuringStimulus" and "AfterStimulus" with corresponding values
indicating valid arguments for matplotlib color.
"mplPars": dict, containing matplotlib rcParams to be overridden.

Usage: python makeContImages <json param file>
'''

import sys
import os
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import quantities as qu
from GJEphys.rawDataAnalyse import RawDataAnalyser
import json

assert len(sys.argv) == 2, 'Improper Usage! Please use as follows:\npython makeContImages <json param file>'

with open(sys.argv[1]) as fle:

    pars = json.load(fle)
    NIXPath = pars['NIXPath']
    expName = pars['expName']
    freqs = pars['freqs']
    catResDir = pars['catResDir']
    downSampleFactor = pars['downSampleFactor']
    type2Color = pars['type2Color']
    mplPars = pars['mplPars']

plt.rcParams.update(mplPars)

plt.ioff()

rda = RawDataAnalyser(dirpath=NIXPath, expName=expName)

resps = rda.getContResps(freqs=freqs)
spikes = rda.getContSpikes(freqs=freqs)

expDir = os.path.join(catResDir, expName)
if not os.path.exists(expDir):
    os.mkdir(expDir)

for freq, freqResps in resps.iteritems():

    print('Doing {}'.format(freq))

    overallFig, overallAx = plt.subplots(nrows=1, ncols=1, figsize=(7, 5.6))

    freqDir = os.path.join(expDir, str(freq))
    if not os.path.exists(freqDir):
        os.mkdir(freqDir)

    nTrial = len(freqResps)
    trialStarts = []
    bottoms = [0 * qu.mV for x in range(nTrial + 1)]
    zeros = []

    for trialInd, trialResps in enumerate(freqResps):

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5.6))
        fig.suptitle('{} Hz; Trial no:{}'.format(freq, trialInd + 1))

        trialMax = reduce(max, [x.max() for x in trialResps.itervalues()])
        trialMin = reduce(min, [x.min() for x in trialResps.itervalues()])


        trialStart = trialResps['DuringStimulus'].t_start
        trialStart.units = qu.s

        trialStarts.append(trialStart.magnitude)

        offset = bottoms[trialInd] - trialMin

        zeros.append(offset.magnitude)

        for typ, typeResp in trialResps.iteritems():
            
            ax[0].plot(typeResp.times[::downSampleFactor], typeResp[::downSampleFactor], color=type2Color[typ],
                               marker=None, ls='-')

            overallAx.plot(typeResp.times[::downSampleFactor] - trialStart,
                           typeResp[::downSampleFactor] + offset,
                           color=type2Color[typ], marker=None, ls='-')


            typeSpikes = spikes[freq][trialInd][typ]
            typeSpikesTimes = qu.Quantity(typeSpikes.times)
            ax[0].plot(typeSpikesTimes, [max(typeResp)] * len(typeSpikes), 'k^')

            overallAx.plot(typeSpikesTimes - trialStart, [max(typeResp) + offset] * len(typeSpikes), 'k^')

            if typ == 'DuringStimulus':
                ax[1].plot(typeResp.times[::downSampleFactor], typeResp[::downSampleFactor], color=type2Color[typ],
                           marker=None, ls='-')
                ax[1].plot(typeSpikesTimes, [max(typeResp)] * len(typeSpikes), 'k^')

        bottoms[trialInd + 1] = bottoms[trialInd] + 1.1 * (trialMax - trialMin)

        ax[1].set_xlabel('time (s)')
        fig.text(0, 0.5, 'Membrane Potential (mV)', va='center', rotation='vertical')
        fig.tight_layout()
        fig.canvas.draw()

        fig.savefig(os.path.join(freqDir, 'Trial{}.png'.format(trialInd)), dpi=300)

        plt.close(fig.number)
        del fig


    overallAx.set_xlabel('time relative to \n application of stimulus (s)')
    overallAx.set_ylabel('time of application \n of stimulus (s)')
    overallAx.set_yticks(zeros)
    overallAx.set_yticklabels(map(str, trialStarts))
    overallFig.tight_layout()
    overallFig.savefig(os.path.join(expDir, 'overall{:d}Hz.png'.format(int(freq))), dpi=300)

    plt.close(overallFig.number)
    del overallFig

del resps
del spikes
del rda
