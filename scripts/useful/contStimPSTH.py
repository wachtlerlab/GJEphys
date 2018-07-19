'''
This script has functions to save the spike data of all trials of a set of specified experiments and plot a PSTH of
the responses of those trials from the saved data.
The script calls different function depending on the first command line argument.
Please look at the documentation of individual functions for more information.

Usage: Execute this script without any command line arguments to get a list of usages
'''

import pandas as pd
import quantities as qu
from GJEphys.rawDataAnalyse import RawDataAnalyser
import numpy as np
from matplotlib import pyplot as plt
import sys
from GJEphys.pdColumnNameMapCont import mdFN
from GJEphys.matplotlibRCParams import mplPars
import seaborn as sns
from GJEphys.ssvkernel import ssvkernel
import warnings
warnings.filterwarnings(action='once')


def saveContStimPSTH(inputXL, dataXL):
    '''
    Extracts the spike times of each trial of the set of experiments specified in <inputXL> and saves them into the
    excel file <dataXL> with one spike data per row, including the metadata like <Experiment ID>, <Labor State>
    and <TrialName>
    :param inputXL: string, path of an excel file. The excel file should contain three columns with headings
    "Experiment ID", "Labor State" and "NIX File Directory".
    :param dataXL: string, path where resulting excel file will be written.
    :return:
    '''
    types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']
    typeDurs = [3 * qu.s, 1 * qu.s, 1 * qu.s]

    inputDF = pd.read_excel(inputXL)

    dataDF = pd.DataFrame()

    for ind, (expID, laborState, nixPath) in inputDF.loc[:,
                                             ["Experiment ID", "Labor State", "NIX File Directory"]].iterrows():
        rda = RawDataAnalyser(expID, nixPath)

        expSpikes = rda.getContSpikes(types=types, freqs=[265])

        for trialInd, trialSpikes in enumerate(expSpikes[265]):
            print("Doing {}, Trial {}".format(expID, trialInd + 1))
            trialSpikeTimeList = []
            for typeSpikeTrain in trialSpikes.itervalues():
                temp = typeSpikeTrain.copy()
                temp.units = qu.s
                trialSpikeTimeList.append(temp.times.magnitude)
            allTrialSpikeTimes = np.concatenate(trialSpikeTimeList)

            spikeTimeUnits = qu.s
            allTrialSpikeTimes *= spikeTimeUnits

            allTrialSpikeTimes -= trialSpikes["DuringStimulus"].t_start

            allTrialSpikeTimes = allTrialSpikeTimes[
                np.logical_and(-typeDurs[0] <= allTrialSpikeTimes, allTrialSpikeTimes <= typeDurs[1] + typeDurs[2])]

            for spTime in allTrialSpikeTimes:
                tempS = pd.Series()
                tempS[mdFN["expID"]] = expID
                tempS[mdFN["freq"]] = 265
                tempS[mdFN["laborState"]] = laborState
                tempS[mdFN["trialName"]] = "Trial{}".format(trialInd)
                tempS[mdFN["trialStart"]] = trialSpikes["DuringStimulus"].t_start
                tempS["Spike Time (s)"] = spTime

                dataDF = dataDF.append(tempS, ignore_index=True)

    dataDF.to_excel(dataXL)


def saveFRTS(dataXL, statsXL):
    '''
    Uses spike data stored in <dataXL> to calculate firing rate for each trial and for time bins 20ms long over the
    period starting 500ms before stimulus onset and ending 500ms after stimulus offset. Resulting values are stored one
    firing rate per row along with metadata like <time (s)> (center of time bin), <Trial Name> and <Experiment ID>.
    :param dataXL: string, path to an excel file, generated using the function saveData above
    :param statsXL: string, path where the result data will be stored.
    :return:
    '''
    dataDF = pd.read_excel(dataXL)


    statsDF = pd.DataFrame()

    binEdges = np.arange(-0.5, 1.5, 0.02)
    binWidth = binEdges[1] - binEdges[0]
    binCenters = binEdges[:-1] + 0.5 * binWidth

    for expID, expIDDF in dataDF.groupby(mdFN["expID"]):
        for trialInd, trialDF in expIDDF.groupby(expIDDF[mdFN["trialName"]]):
            print("Doing {} {}".format(expID, trialInd))
            counts, bins = np.histogram(trialDF["Spike Time (s)"], binEdges)

            for count, binC in zip(counts, binCenters):

                firingRate = count / binWidth
                tempS = pd.Series()
                tempS["Firing Rate (spikes/s)"] = firingRate
                tempS["time (s)"] = binC
                tempS[mdFN["expID"]] = expID
                tempS[mdFN["trialName"]] = trialInd
                lsUnique = trialDF[mdFN["laborState"]].unique()
                assert len(lsUnique) == 1, "Inconsistent labor state entries for {}, {}".format(expID, trialInd)
                tempS[mdFN["laborState"]] = lsUnique[0]
                statsDF = statsDF.append(tempS, ignore_index=True)


    statsDF.to_excel(statsXL)


def plotFRTS(statsXL, outBase):
    """
    Plots average firing rate vs time bar plot from <statsXL>, separating and averaging over each unique value in the
    column <Labor State>.
    :param statsXL: string, path to an excel file, generated using "saveFRTS" function
    :param outBase: string, figure is saved as "<outBase>.png"
    :return:
    """
    statsDF = pd.read_excel(statsXL)

    binEdges = np.arange(-0.5, 1.5, 0.02)
    binWidth = binEdges[1] - binEdges[0]
    binCenters = binEdges[:-1] + 0.5 * binWidth

    mplPars["lines.linewidth"] = 0.7
    sns.set(rc=mplPars, style='darkgrid')
    fig1, ax1 = plt.subplots(figsize=(7, 5.6))

    if len(statsDF[mdFN["laborState"]].unique()) > 1:
        sns.barplot(data=statsDF, x="time (s)", y="Firing Rate (spikes/s)", hue=mdFN["laborState"], ax=ax1, dodge=True,
                  ci=None, palette=["b", "r"])
    else:
        sns.barplot(data=statsDF, x="time (s)", y="Firing Rate (spikes/s)", ax=ax1,
                      ci=None, color='b')
    currentXTicks = np.round(binCenters, 2).tolist()
    ax1.set_xticks([currentXTicks.index(0.01) - 0.5, currentXTicks.index(1.01) - 0.5])
    ax1.set_xticklabels(['Stimulus Onset', 'Stimulus Offset'])
    ax1.set_ylabel("Average Firing\nRate (spikes/s)")
    ax1.set_xlabel("")
    ax1.grid(True)
    fig1.tight_layout()
    fig1.savefig("{}.png".format(outBase), dpi=300)


def plotFRTSKDE(dataXL, outBase):
    """
    Estimates the firing rate profile from the spiking data in <dataXL> using AdaptiveKDE[1][2], plots and saves it.
    :param dataXL: string, path of the excel file generated by the function "saveContStimPSTH".
    :param outBase: string, output file will be saved as "<outBase>.png"
    :return:

    References
    [1] Shimazaki & Shinomoto "Kernel bandwidth optimization in spike rate estimation" J Comput Neurosci (2010) 29: 171.
    https://doi.org/10.1007/s10827-009-0180-4
    [2] https://github.com/cooperlab/AdaptiveKDE
    """

    dataDF = pd.read_excel(dataXL)
    plt.rcParams.update(mplPars)

    binWidth = 0.001
    tEdges = np.arange(-1, 2, binWidth)
    tEst = tEdges[1:] - 0.5 * binWidth

    mplPars["lines.linewidth"] = 0.7
    sns.set(rc=mplPars, style='ticks')
    fig1, ax1 = plt.subplots(figsize=(7, 5.6))
    fig2, ax2 = plt.subplots(figsize=(7, 5.6))
    fig3, ax3 = plt.subplots(figsize=(7, 5.6))

    nLS = len(dataDF[mdFN["laborState"]].unique())

    if nLS == 1:
        palette = ["b"]
    elif nLS == 2:
        palette = ["b", "r"]
    else:
        palette = sns.color_palette("RdYlBu")
    for lsInd,  (ls, lsDF) in enumerate(dataDF.groupby(mdFN["laborState"])):
        spikeData = lsDF
        spikeTimes = spikeData["Spike Time (s)"]
        trialExp = spikeData.apply(lambda x: "{}{}".format(x[mdFN["expID"]], x[mdFN["trialName"]]), axis=1)
        nTrials = trialExp.unique().size
        nSpikes = spikeData.shape[0]
        res = ssvkernel(x=spikeTimes.values, tin=tEst, WinFunc="Gauss")

        fr = res[0] * nSpikes / nTrials
        confInt95 = res[5][0] * nSpikes / nTrials
        optW = res[2]

        ax1.plot(tEst, fr, color=palette[lsInd], ls="-", marker="None", label=ls)
        ax1.fill_between(tEst, fr - 0.5 * confInt95, fr + 0.5 * confInt95, facecolors=palette[lsInd], interpolate=True,
                         alpha=0.2)

        ax2.plot(tEst, fr, color=palette[lsInd], ls="-", marker="None", label=ls)
        ax2.fill_between(tEst, fr - 0.5 * confInt95, fr + 0.5 * confInt95, facecolors=palette[lsInd], interpolate=True,
                         alpha=0.2)

        ax3.plot(tEst, optW, color=palette[lsInd], ls="-", marker="None", label=ls)

    ax1.set_ylabel("Average Firing\nRate(spikes/s)")
    ax1.set_xlabel("time (s)")

    ax1.set_xticks(np.arange(-1, 2, 0.1), minor=True)
    ax1.grid(b=True, axis='x', which="major", linestyle='--')
    ax1.grid(b=False, axis='y')
    ax1.tick_params(which="both", direction="in", bottom=True, top=True)

    ax2.set_ylabel("Average Firing\nRate(spikes/s)")
    ax2.set_xlabel("time (s)")
    ax2.set_xticks(np.arange(-1, 2, 0.01), minor=True)
    ax2.grid(b=True, axis='x', which="major", linestyle='--')
    ax2.grid(b=False, axis='y')
    ax2.tick_params(which="both", direction="in", bottom=True, top=True)

    ax3.set_ylabel("Std. Dev. of\nSmoothing Gaussian (s)")
    ax3.set_xlabel("time (s)")
    ax3.tick_params(direction="in")



    if nLS > 1:
        ax1.legend(loc="upper center", ncol=2)
        # ax2.legend(loc="upper center", ncol=2)

    fig1.tight_layout()
    fig1.savefig("{}.png".format(outBase), dpi=300)

    ax2.set_xlim((-0.05, 0.1))
    fig2.tight_layout()
    fig2.savefig("{}_OnPhasic.png".format(outBase), dpi=300)

    fig3.tight_layout()
    fig3.savefig("{}_optW.png".format(outBase), dpi=300)


def plotContStimPSTH(dataXL, outBase):
    """
    Plots a PSTH plot, combining data from all trials in <dataXL>.
    :param dataXL: string, path of an excel file.
    :param outBase: string, figure is saved as "<outBase>.png"
    :return:
    """
    dataDF = pd.read_excel(dataXL)
    plt.rcParams.update(mplPars)

    for expID, expIDDF in dataDF.groupby(mdFN["expID"]):
        print("Doing {}".format(expID))
        fig, ax = plt.subplots(figsize=(7, 5.6))
        ax.hist(expIDDF["Spike Time (s)"], bins=np.arange(-0.5, 1.5, 0.01), density=True, stacked=True,
                 histtype="stepfilled")
        ax.set_xlabel("time relative to stimulus onset(s)")
        ax.set_ylabel("Normalized Frequency")
        nTrials = expIDDF[mdFN["trialName"]].unique().shape[0]
        ax.set_title("Continuous stimulus PSTH\n for {} (nTrials={})".format(expID, nTrials))
        xticks = np.arange(-0.45, 1.5, 0.05)
        ax.set_xticks(xticks)
        ax.set_xticklabels(["{:1.2f}".format(xticks[x]) if x % 2 else '' for x in range(xticks.shape[0])], rotation=90)
        ax.grid(True, axis='x')
        fig.tight_layout()
        fig.savefig("{}_{}.png".format(outBase, expID), dpi=300)

    fig1, ax1 = plt.subplots(figsize=(7, 5.6))
    ax1.hist(dataDF["Spike Time (s)"], bins=np.arange(-0.5, 1.5, 0.01), density=True, stacked=True,
             histtype="stepfilled")
    ax1.set_xlabel("time relative to stimulus onset(s)")
    ax1.set_ylabel("Normalized Frequency")
    ax1.set_title("Continuous stimulus PSTH for DL-Int-1")
    xticks = np.arange(-0.45, 1.5, 0.05)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(["{:1.2f}".format(xticks[x]) if x % 2 else '' for x in range(xticks.shape[0])], rotation=90)
    ax1.grid(True, axis='x')
    fig1.tight_layout()
    fig1.savefig("{}.png".format(outBase), dpi=300)


if __name__ == "__main__":

    assert len(sys.argv) == 4, "Improper Usage! Please use as\n" \
                               "python {currFile} save <input excel file> " \
                               "<data excel file> or \n" \
                               "python {currFile} plot <data excel file> " \
                               "<output base> or \n" \
                               "python {currFile} saveFRTS <data excel file> " \
                               "<output statsXL> or \n" \
                               "python {currFile} plotFRTS <data excel file> " \
                               "<output base>or \n" \
                               "python {currFile} plotFRTSKDE <data excel file> " \
                               "<output base>".format(currFile=sys.argv[0])

    if sys.argv[1] == "save":
        saveContStimPSTH(*sys.argv[2:])
    elif sys.argv[1] == "plot":
        plotContStimPSTH(*sys.argv[2:])
    elif sys.argv[1] == "saveFRTS":
        saveFRTS(*sys.argv[2:])
    elif sys.argv[1] == "plotFRTS":
        plotFRTS(*sys.argv[2:])
    elif sys.argv[1] == "plotFRTSKDE":
        plotFRTSKDE(*sys.argv[2:])

