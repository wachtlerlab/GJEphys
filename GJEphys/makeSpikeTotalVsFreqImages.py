from GJEphys.rawDataAnalyse import RawDataAnalyser
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import quantities as qu
import os
import sys
import json
import numpy as np

def makeSpikeTotalVsFreqImages(expNamesDict, NIXPath, freqs, resDir, mplPars):

    sns.set(style="whitegrid", rc=mplPars)

    assert len(expNamesDict) == 1, "expNamesDict must be a dict of one key-value pair!"

    categoryStatsDF = pd.DataFrame()
    for expName in expNamesDict.values()[0]:

        print("Doing {}".format(expName))
        rda = RawDataAnalyser(dirpath=NIXPath, expName=expName)
        spikes = rda.getContSpikes(freqs=freqs)

        expStatsDF = pd.DataFrame()

        for freq, freqSpikes in spikes.iteritems():

            tempS = pd.Series()

            for trialInd, trialSpikesDict in enumerate(freqSpikes):

                trialSpikes = trialSpikesDict["DuringStimulus"]
                tempS["Experiment ID"] = expName
                tempS["Frequency (Hz)"] = freq
                tempS["Trial Number"] = trialInd
                tempS["Neuron Category"] = expNamesDict.keys()[0]
                stimDur = trialSpikes.duration
                stimDur.units = qu.s
                tempS["Spike Rate (spikes per second)"] = trialSpikes.shape[0] / stimDur.magnitude

                expStatsDF = expStatsDF.append(tempS, ignore_index=True)

        categoryStatsDF = categoryStatsDF.append(expStatsDF, ignore_index=True)

    fig, ax = plt.subplots(figsize=(7, 5.6))
    fg = sns.factorplot(data=categoryStatsDF, x="Frequency (Hz)", y="Spike Rate (spikes per second)", hue="Experiment ID",
                   ax=ax, kind='point', ci=95, n_boot=1000)
    ax.legend(bbox_to_anchor=(1.65, 1))

    ax.set_title(expNamesDict.keys()[0])
    fig.tight_layout(rect=[0, 0, 0.75, 0.9])
    fig.savefig(os.path.join(resDir, "IndividualExperimentsSeparately.png"), dpi=300)

    fig1, ax1 = plt.subplots(figsize=(7, 5.6))
    sns.violinplot(data=categoryStatsDF, x="Frequency (Hz)", y="Spike Rate (spikes per second)",
                   ax=ax1, scale="area", scale_hue=False)
    ax1.set_title(expNamesDict.keys()[0])
    fig1.tight_layout()
    fig1.savefig(os.path.join(resDir, "AllExperimentsCombined.png"), dpi=300)


if __name__ == "__main__":

    assert len(sys.argv) == 2, 'Improper Usage! Please use as follows:\npython makeContImages <json param file>'

    with open(sys.argv[1]) as fle:
        pars = json.load(fle)
        expNamesDict = pars['expNamesDict']
        NIXPath = pars["NIXPath"]
        freqs = pars['freqs']
        catResDir = pars['catResDir']
        mplPars = pars['mplPars']

    makeSpikeTotalVsFreqImages(expNamesDict=expNamesDict,
                               NIXPath=NIXPath,
                               freqs=freqs,
                               resDir=catResDir,
                               mplPars=mplPars)




