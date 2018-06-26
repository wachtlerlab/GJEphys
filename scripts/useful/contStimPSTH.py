import pandas as pd
import quantities as qu
from GJEphys.rawDataAnalyse import RawDataAnalyser
import numpy as np
from matplotlib import pyplot as plt
import sys
from GJEphys.pdColumnNameMapCont import mdFN
from GJEphys.matplotlibRCParams import mplPars
import seaborn as sns



def saveContStimPSTH(inputXL, dataXL):
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

    statsDF = pd.read_excel(statsXL)

    binEdges = np.arange(-0.5, 1.5, 0.02)
    binWidth = binEdges[1] - binEdges[0]
    binCenters = binEdges[:-1] + 0.5 * binWidth

    mplPars["lines.linewidth"] = 0.7
    sns.set(rc=mplPars, style='darkgrid')
    fig1, ax1 = plt.subplots(figsize=(7, 5.6))

    if len(statsDF[mdFN["laborState"]].unique()) > 1:
        sns.barplot(data=statsDF, x="time (s)", y="Firing Rate (spikes/s)", hue=mdFN["laborState"], ax=ax1, dodge=True,
                  ci=None)
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


def plotContStimPSTH(dataXL, outBase):

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
                               "<output base>".format(currFile=sys.argv[0])

    if sys.argv[1] == "save":
        saveContStimPSTH(*sys.argv[2:])
    elif sys.argv[1] == "plot":
        plotContStimPSTH(*sys.argv[2:])
    elif sys.argv[1] == "saveFRTS":
        saveFRTS(*sys.argv[2:])
    elif sys.argv[1] == "plotFRTS":
        plotFRTS(*sys.argv[2:])

