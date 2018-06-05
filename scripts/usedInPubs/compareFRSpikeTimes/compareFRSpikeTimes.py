'''
Description: This script contains functions for extracting and plotting the following features of DL-Int-1 ephys response:
1. firing rates in four different phases of DL-Int-1 ephys reponse
2. spike timing of the first four spikes after stimulus onset
The script calls different function depending on the first command line argument.
Please look at the documentation of individual functions for more information.

Usage: Execute this script without any command line arguments to get a list of usages.
'''


from GJEphys.pdColumnNameMapCont import mdFN, spikeFRSpikeTimesFNs, spikeFRSpikeTimesFuncs
from GJEphys.rawDataAnalyse import RawDataAnalyser
import pandas as pd
import quantities as qu
import numpy as np
from neo import SpikeTrain
import sys
from GJEphys import additionalEphysFuncs
import matplotlib.pyplot as plt
import seaborn as sns
from GJEMS.viz.matplotlibRCParams import mplPars
from scipy.stats import ttest_ind
import os
import shutil

def saveData(inputXL, dataXL):
    '''
    Extracts the following two features from the ephys responses of DL-Int-1:
     1. firing rates in four different phases of DL-Int-1 ephys reponse
     2. spike timing of the first four spikes after stimulus onset
    Features are extracted for the Experiment IDs specified in 'inputXL'. For each row of 'inputXL', a NIX file
    with the name "<Experiment ID>.h5" in the directory <NIX File Directory> is used to extract features.
    The extracted feature values are stored in the files specified by 'dataXL'.
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

            trialSpikeTrain = SpikeTrain(times=allTrialSpikeTimes,
                                         units=spikeTimeUnits,
                                         t_start=-typeDurs[0], t_stop=typeDurs[1] + typeDurs[2])
            tempS = pd.Series()
            tempS[mdFN["expID"]] = expID
            tempS[mdFN["freq"]] = 265
            tempS[mdFN["laborState"]] = laborState
            tempS[mdFN["trialName"]] = "Trial{}".format(trialInd)
            tempS[mdFN["trialStart"]] = trialSpikes["DuringStimulus"].t_start

            for funcName, funcNameFull in spikeFRSpikeTimesFNs.iteritems():

                func = getattr(additionalEphysFuncs, spikeFRSpikeTimesFuncs[funcName])
                tempS[funcNameFull] = func(None, trialSpikeTrain, None, None)

            dataDF = dataDF.append(tempS, ignore_index=True)

    dataDF.to_excel(dataXL)


def plotFRData(dataXL, outBase):
    '''
    Uses firing rate feature values extracted in the function 'saveData' above and makes a violin plot using seaborn
    with the four response interval types on the X-Axes and Firing rate on Y-Axis.
    The figure generated is saved to the file "<outBase>_FR.png"
    :param dataXL: string, path of the excel file generated using the function 'saveData' above.
    :param outBase: string, the generated plot will be saved to the file "<outBase>_FR.png"
    :return:
    '''

    dataDF = pd.read_excel(dataXL)

    FRColKeys = ["spontFR3", "initFR", "laterFR", "reboundFR"]
    FRCols = [spikeFRSpikeTimesFNs[x] for x in FRColKeys]

    FRData = dataDF.loc[:, FRCols+[mdFN["laborState"]]]
    FRData.set_index(mdFN["laborState"], inplace=True)
    FRData.rename(columns={spikeFRSpikeTimesFNs[k]: k for k in FRColKeys}, inplace=True)
    FRDataStacked = FRData.stack().reset_index()
    FRDataStacked.rename(columns={0: "Firing Rate (spikes/s)",
                                  "level_1": "Interval"}, inplace=True)

    sns.set(style="darkgrid", rc=mplPars)

    fig, ax = plt.subplots(figsize=(7, 5.6))

    sns.violinplot(hue=mdFN['laborState'], y="Firing Rate (spikes/s)", x="Interval",
                data=FRDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=["r", "b"],
                order=FRColKeys, split=True, scale="area", inner=None, cut=0, linewidth=0, scale_hue=False)

    sns.pointplot(hue=mdFN['laborState'], y="Firing Rate (spikes/s)", x="Interval",
                data=FRDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=['w', 'w'], join=False,
                  markers="_", dodge=True, ci=None, size=30)

    for colInd, (col, s) in enumerate(FRData.iteritems()):
        t, pVal = ttest_ind(s["Forager"], s["Newly Emerged"], equal_var=False)


        col = 'k'
        if pVal < 0.05:
            col = 'g'

        ax.text(colInd, -5, "{:1.1e}".format(pVal), fontdict={'color': col}, fontsize=plt.rcParams['xtick.labelsize'],
                horizontalalignment='center', verticalalignment='center')


    ax.set_xticklabels(["Spontaneous\nActivity (3s)", "On-phasic\nResponse (75ms)",
                        "Sustained\nResponse (925ms)", "Off-phasic\nResponse (50ms)"], rotation=45)
    ax.set_ylim(-10, 80)
    ax.set_xlabel("")

    l1, = ax.plot((), (), 'rs', label="Newly Emerged")
    l2, = ax.plot((), (), 'bs', label="Forager")
    ax.legend(handles=[l1, l2], loc="upper right")

    fig.tight_layout()
    fig.savefig("{}_FR.png".format(outBase), dpi=150)


def plotSpikeTimes(dataXL, outBase):
    '''
    Uses the spike time feature values extracted by 'saveData' above and makes a violin plot using seaborn with the
    type of spike timing on X-Axis and spike time values on Y-axis.
    :param dataXL: string, path of the excel file generated using the function 'saveData' above.
    :param outBase: string, the generated plot will be saved to the file "<outBase>_spikeTimes.png"
    :return:
    '''

    dataDF = pd.read_excel(dataXL)

    spikeTimesColKeys = ["firstSpikeLatency", "secondSpikeBISI",
                         "thirdSpikeBISI", "fourthSpikeBISI"]
    spikeTimesCols = [spikeFRSpikeTimesFNs[x] for x in spikeTimesColKeys]

    spikeTimesData = dataDF.loc[:, spikeTimesCols+[mdFN["laborState"]]]
    spikeTimesData.set_index(mdFN["laborState"], inplace=True)
    spikeTimesDataStacked = spikeTimesData.stack().reset_index()
    spikeTimesDataStacked.rename(columns={0: "Time (ms)",
                                  "level_1": "Interval"}, inplace=True)

    sns.set(style="darkgrid", rc=mplPars)

    fig, ax = plt.subplots(figsize=(7, 5.6))

    sns.violinplot(hue=mdFN['laborState'], y="Time (ms)", x="Interval",
                data=spikeTimesDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=["r", "b"],
                order=spikeTimesCols, split=True, scale="area", inner=None, cut=0, linewidth=0, scale_hue=False)

    sns.pointplot(hue=mdFN['laborState'], y="Time (ms)", x="Interval",
                data=spikeTimesDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=['w', 'w'], join=False,
                  markers="_", dodge=True, ci=None, size=30)

    for colInd, (col, s) in enumerate(spikeTimesData.iteritems()):
        t, pVal = ttest_ind(s["Forager"], s["Newly Emerged"], equal_var=False, nan_policy='omit')
        col = 'k'
        if pVal < 0.05:
            col = 'g'

        ax.text(colInd, -5, "{:1.1e}".format(pVal), fontdict={'color': col}, fontsize=plt.rcParams['xtick.labelsize'],
                horizontalalignment='center', verticalalignment='center')

    ax.set_xticklabels(spikeTimesCols, rotation=45)
    ax.set_ylim(-10, 80)
    ax.set_xlabel("")

    l1, = ax.plot((), (), 'rs', label="Newly Emerged")
    l2, = ax.plot((), (), 'bs', label="Forager")
    ax.legend(handles=[l1, l2], loc="upper right")

    fig.tight_layout()
    fig.savefig("{}_spikeTimes.png".format(outBase), dpi=150)


def plotInhFRVsSpontFR(dataXL, outDir):
    """
    Creates a directory at <outDir> is none exists already. Uses firing rate feature values extracted in the function
    'saveData' above. For every unique <Experiment ID>, plots firing rate values of spontaneous activity
    in an interval 3s before stimulus onset on X axis and firing rate values during sustained response (75-1000ms) on
    Y-axis. The resulting figures are saved in <outDir> under the name '<expID>.png'.  A few other plots collecting
    such points from all <Experiment ID>s and using different coloring schemes are also generated in <outDir>.
    :param dataXL: string, path of the excel file generated using the function 'saveData' above.
    :param outDir: string, directory where generated figures are saved as images
    :return:
    """

    sns.set(style="darkgrid", rc=mplPars)

    dataDF = pd.read_excel(dataXL)

    if os.path.isdir(outDir):
        shutil.rmtree(outDir)
    os.mkdir(outDir)

    nExpIDs = len(dataDF["Experiment ID"].unique())

    allFig, allAx = plt.subplots(figsize=(7, 5.6))
    fvneFig, fvneAx = plt.subplots(figsize=(7, 5.6))


    for expInd, (expID, expDF) in enumerate(dataDF.groupby("Experiment ID")):

        print("Doing {}".format(expID))

        expCol = plt.cm.spectral((expInd + 1) / float(nExpIDs))



        allAx.plot(expDF[spikeFRSpikeTimesFNs["spontFR3"]],
                   expDF[spikeFRSpikeTimesFNs["laterFR"]], color=expCol, ms=5, marker="o", ls="None")

        if expDF.size:
            if expDF.iloc[0, :][mdFN["laborState"]] == "Forager":
                lsCol = "b"
            else:
                lsCol = 'r'
            fvneAx.plot(expDF[spikeFRSpikeTimesFNs["spontFR3"]],
                       expDF[spikeFRSpikeTimesFNs["laterFR"]], color=lsCol, ms=5, marker="o", ls="None", label=None)

        fig, ax = plt.subplots(figsize=(7, 5.6))
        maxTrials = max(float(s[5:]) + 1 for s in expDF["TrialName"])
        for rowID, rowS in expDF.iterrows():
            trialNo = float(rowS["TrialName"][5:]) + 1
            col = plt.cm.spectral(1 - trialNo / maxTrials)
            ax.plot(rowS[spikeFRSpikeTimesFNs["spontFR3"]],
                    rowS[spikeFRSpikeTimesFNs["laterFR"]],
                    color=col, marker="o", ms=10)

        ax.set_ylabel("FR during Sustained \n Response (75 to 1000ms)")
        ax.set_xlabel("FR during Spontaneous \n Response (-3000 to 0ms)")
        ax.set_title("VIBGYOR colormap\nViolet=Trial {}, Red=Trial {}".format(1, maxTrials))

        fig.tight_layout()
        fig.savefig(os.path.join(outDir, "{}.png".format(expID)), dpi=300)
        plt.close(fig.number)

    allAx.set_ylabel("FR during Sustained \n Response (75 to 1000ms)")
    allAx.set_xlabel("FR during Spontaneous \n Response (-3000 to 0ms)")
    maxX = 25
    allAx.plot((0, maxX), (0, maxX), 'k:')
    allAx.set_ylim(-5, 60)
    allAx.set_xlim(-5, 30)

    allFig.tight_layout()
    allFig.savefig(os.path.join(outDir, "allExps.png"), dpi=300)

    allAx.set_xlim(-5, 15)
    allAx.set_ylim(-5, 20)

    allFig.tight_layout()
    allFig.savefig(os.path.join(outDir, "allExps_zoomed.png"), dpi=300)

    fvneAx.set_ylabel("FR during Sustained \n Response (75 to 1000ms)")
    fvneAx.set_xlabel("FR during Spontaneous \n Response (-3000 to 0ms)")

    fvneAx.plot((0, maxX), (0, maxX / 3), 'k:')
    fvneAx.set_ylim(-5, 60)
    fvneAx.set_xlim(-5, 30)

    l1, = fvneAx.plot((), (), 'bo', ms=5)
    l2, = fvneAx.plot((), (), 'ro', ms=5)
    fvneAx.legend((l1, l2), ("Forager", "Newly Emerged"), loc="best")
    fvneFig.tight_layout()
    fvneFig.savefig(os.path.join(outDir, "allExpsFVsNE.png"), dpi=300)


def plotFRRelative2SpontFR(dataXL, outBase):
    """
    Generates images with figures similar to those generated by the function "plotFR" above. From all firing rate
    values, corresponding spontaneous activity rates are subtracted before plotting.
    :param dataXL: string, path of the excel file generated using the function 'saveData' above.
    :param outBase: string, the generated plot will be saved to the file "<outBase>_spikeTimes.png"
    :return:
    """

    sns.set(style="darkgrid", rc=mplPars)

    dataDF = pd.read_excel(dataXL)

    dataDF["initFRRel2spontFR"] = dataDF[spikeFRSpikeTimesFNs["initFR"]] - dataDF[spikeFRSpikeTimesFNs["spontFR3"]]
    dataDF["laterFRRel2spontFR"] = dataDF[spikeFRSpikeTimesFNs["laterFR"]] - dataDF[spikeFRSpikeTimesFNs["spontFR3"]]
    dataDF["reboundFRRel2spontFR"] = dataDF[spikeFRSpikeTimesFNs["reboundFR"]] \
                                     - dataDF[spikeFRSpikeTimesFNs["spontFR3"]]

    FRColKeys = ["initFRRel2spontFR", "laterFRRel2spontFR", "reboundFRRel2spontFR"]

    FRData = dataDF.loc[:, FRColKeys + [mdFN["laborState"]]]
    FRData.set_index(mdFN["laborState"], inplace=True)
    FRDataStacked = FRData.stack().reset_index()
    FRDataStacked.rename(columns={0: "Firing Rate relative to \n"
                                     "Spontaneous activity\n(spikes/s)",
                                  "level_1": "Interval"}, inplace=True)

    sns.set(style="darkgrid", rc=mplPars)

    fig, ax = plt.subplots(figsize=(7, 5.6))

    sns.violinplot(hue=mdFN['laborState'], y="Firing Rate relative to \nSpontaneous activity\n(spikes/s)", x="Interval",
                   data=FRDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=["r", "b"],
                   order=FRColKeys, split=True, scale="area", inner=None, cut=0, linewidth=0, scale_hue=False)

    sns.pointplot(hue=mdFN['laborState'], y="Firing Rate relative to \nSpontaneous activity\n(spikes/s)", x="Interval",
                  data=FRDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=['w', 'w'], join=False,
                  markers="_", dodge=True, ci=None, size=30)

    for colInd, (col, s) in enumerate(FRData.iteritems()):
        t, pVal = ttest_ind(s["Forager"], s["Newly Emerged"], equal_var=False)

        col = 'k'
        if pVal < 0.05:
            col = 'g'

        ax.text(colInd, -5, "{:1.1e}".format(pVal), fontdict={'color': col}, fontsize=plt.rcParams['xtick.labelsize'],
                horizontalalignment='center', verticalalignment='center')

    ax.set_xticklabels(["On-phasic\nResponse (75ms)",
                        "Sustained\nResponse (925ms)", "Off-phasic\nResponse (50ms)"], rotation=45)
    ax.set_ylim(-40, 80)
    ax.set_xlabel("")

    l1, = ax.plot((), (), 'rs', label="Newly Emerged")
    l2, = ax.plot((), (), 'bs', label="Forager")
    ax.legend(handles=[l1, l2], loc="upper right")

    fig.tight_layout()
    fig.savefig("{}_FRRelative2spontFR.png".format(outBase), dpi=150)



if __name__ == "__main__":

    assert len(sys.argv) == 4, "Improper Usage! Please use as\n" \
                               "python {currFile} save <input excel file> " \
                               "<output data file>  or \n" \
                               "python {currFile} plotFR <input Data excel file> " \
                               "<output Base> or\n" \
                               "python {currFile} plotSpikeTimes <input Data excel file> " \
                               "<output Base> or \n" \
                               "python {currFile} plotLaterFRVsSpontFR <input Data excel file> " \
                               "<output Directory> or \n" \
                               "python {currFile} plotFRRelative2SpontFR <input Data excel file>" \
                               "<output Base>".format(currFile=sys.argv[0])

    if sys.argv[1] == "save":
        saveData(*sys.argv[2:])
    if sys.argv[1] == "plotFR":
        plotFRData(*sys.argv[2:])
    if sys.argv[1] == "plotSpikeTimes":
        plotSpikeTimes(*sys.argv[2:])
    if sys.argv[1] == "plotLaterFRVsSpontFR":
        plotInhFRVsSpontFR(*sys.argv[2:])
    if sys.argv[1] == "plotFRRelative2SpontFR":
        plotFRRelative2SpontFR(*sys.argv[2:])












