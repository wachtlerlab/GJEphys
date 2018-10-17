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
from GJEphys.orphanFuncs import getWelchDF
import pandas as pd
import quantities as qu
import numpy as np
from neo import SpikeTrain
import sys
from GJEphys import additionalEphysFuncs
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from GJEphys.matplotlibRCParams import mplPars, darkTicksPars
from scipy.stats import ttest_ind, linregress, t
import os
import shutil
import warnings
warnings.filterwarnings(action='once')



def saveData(inputXL, dataXL):
    '''
    Extracts the following two features from the ephys responses of DL-Int-1:
     1. firing rates in five different phases of DL-Int-1 ephys reponse
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


def plotFRFVsNE(dataXL, outBase):
    '''
    Uses firing rate feature values extracted in the function 'saveData' above and makes a violin plot using seaborn
    with the four response interval types on the X-Axes and Firing rate on Y-Axis.
    The figure generated is saved to the file "<outBase>_FR.png"
    :param dataXL: string, path of the excel file generated using the function 'saveData' above.
    :param outBase: string, the generated plot will be saved to the file "<outBase>_FR.png"
    :return:
    '''

    dataDF = pd.read_excel(dataXL)

    # FRColKeys = ["spontFR3", "beforeOnsetFR", "initFR", "laterFR", "reboundFR"]
    FRColKeys = ["spontFR3", "initFR", "laterFR", "reboundFR"]
    FRCols = [spikeFRSpikeTimesFNs[x] for x in FRColKeys]

    FRData = dataDF.loc[:, FRCols+[mdFN["laborState"]]]
    FRData.set_index(mdFN["laborState"], inplace=True)
    FRData.rename(columns={spikeFRSpikeTimesFNs[k]: k for k in FRColKeys}, inplace=True)
    FRDataStacked = FRData.stack().reset_index()
    FRDataStacked.rename(columns={0: "Firing Rate (spikes/s)",
                                  "level_1": "Interval"}, inplace=True)

    # hybridPars = dict(mplPars.items() + darkTicksPars.items())
    sns.set(style="whitegrid", rc=mplPars)

    fig, ax = plt.subplots(figsize=(7, 5.6))

    sns.violinplot(hue=mdFN['laborState'], y="Firing Rate (spikes/s)", x="Interval",
                data=FRDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'],
                   palette=[(1, 0, 0, 1), (0, 0, 1, 1)],
                order=FRColKeys, split=True, scale="area", inner=None, cut=0, linewidth=0, scale_hue=False)

    sns.pointplot(hue=mdFN['laborState'], y="Firing Rate (spikes/s)", x="Interval",
                data=FRDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=['w', 'w'], join=False,
                  markers="_", dodge=True, ci=None, size=30)

    for colInd, (col, s) in enumerate(FRData.iteritems()):
        t, pVal = ttest_ind(s["Forager"], s["Newly Emerged"], equal_var=False)


        col = (0, 0, 0, 1)
        if pVal < 0.05:
            col = (0, 0.5, 0, 1)

        ax.text(colInd, -5, "{:1.1e}".format(pVal), fontdict={'color': col}, fontsize=plt.rcParams['xtick.labelsize'],
                horizontalalignment='center', verticalalignment='center')

    # ax.set_xticklabels(["Spontaneous\nActivity (3s)", "Pre-Onset\nActivity (50ms)", "On-phasic\nResponse (75ms)",
    #                     "Inhibitory\nResponse (925ms)", "Rebound\nResponse (75ms)"], rotation=45)
    ax.set_xticklabels(["Spontaneous\nActivity\n(3s)", "On-phasic\nResponse\n(75ms)",
                        "Inhibitory\nResponse\n(925ms)", "Rebound\nResponse\n(75ms)"],
                       rotation=60, va="top", ha="center")
    ax.set_ylim(-10, 80)
    ax.set_xlabel("")

    # l1, = ax.plot((), (), 'rs', label="Newly Emerged")
    # l2, = ax.plot((), (), 'bs', label="Forager")
    # ax.legend(handles=[l1, l2], loc="upper right")

    ax.legend().set_visible(False)

    fig.tight_layout()
    fig.savefig("{}_FR.png".format(outBase), dpi=150)


def plotFRIndNrns(dataXL, outBase):
    '''
    Uses the firing rate values saved by saveData above to plot line plots, one per <Experiment ID>, with firing rates
    on y axis and the five phases on x axis. The figure is saved as "<outBase>.png".
    :param dataXL: string, path to an excel file on file system, generated by saveData above.
    :param outBase: string, used to form the name of the file where output is saved.
    '''
    mplPars["lines.linewidth"] = 0.7
    sns.set(style="darkgrid", rc=mplPars)

    dataDF = pd.read_excel(dataXL)

    FRColKeys = ["spontFR3", "beforeOnsetFR", "initFR", "laterFR", "reboundFR", "afterReboundFR"]
    FRCols = [spikeFRSpikeTimesFNs[x] for x in FRColKeys]

    FRData = dataDF.loc[:, FRCols + [mdFN["expID"]]]
    FRData.set_index(mdFN["expID"], inplace=True)
    FRData.rename(columns={spikeFRSpikeTimesFNs[k]: k for k in FRColKeys}, inplace=True)
    FRDataStacked = FRData.stack().reset_index()
    FRDataStacked.rename(columns={0: "Firing Rate (spikes/s)",
                                  "level_1": "Interval"}, inplace=True)


    expIDLSMap = {}
    for expID, expIDDF in dataDF.groupby(mdFN["expID"]):
        lsUnique = expIDDF[mdFN["laborState"]].unique()
        assert lsUnique.size == 1, "problem with data XL, {} has inconsistent {} entries".format(expID,
                                                                                                 mdFN["laborState"])
        expIDLSMap[expID] = lsUnique[0]

    nForagers = sum([x == "Forager" for x in expIDLSMap.values()])
    nNE = sum([x == "Newly Emerged" for x in expIDLSMap.values()])

    foragerColorMap = sns.light_palette('g', nForagers)
    NEColorMap = sns.light_palette('r', nNE)

    expIDColorMap = {}
    for expID, expIDDF in dataDF.groupby(mdFN["expID"]):
        lsUnique = expIDDF[mdFN["laborState"]].unique()
        assert lsUnique.size==1, "problem with data XL, {} has inconsistent {} entries".format(expID,
                                                                                               mdFN["laborState"])
        if lsUnique[0] == "Forager":
            expIDColorMap[expID] = foragerColorMap.pop(0)
        elif lsUnique[0] == "Newly Emerged":
            expIDColorMap[expID] = NEColorMap.pop(0)



    for expID, expIDDF in FRDataStacked.groupby(mdFN["expID"]):
        fig1, ax1 = plt.subplots(figsize=(7, 5.6))

        sns.pointplot(data=expIDDF, x="Interval", y="Firing Rate (spikes/s)", ci="sd",
                   color=expIDColorMap[expID],
                   ax=ax1)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
        ax1.set_ylim(-10, 100)
        ax1.set_xlabel("")
        ax1.set_title(expID)
        fig1.tight_layout()
        fig1.savefig("{}_indExpIDFR_{}.png".format(outBase, expID), dpi=300)







    fig, ax = plt.subplots(figsize=(7, 5.6))


    sns.pointplot(data=FRDataStacked, x="Interval", y="Firing Rate (spikes/s)", hue=mdFN["expID"], ci="sd",
                   dodge=True, palette=expIDColorMap,
                   ax=ax)

    # for l in ax.lines:
    #     print(l.get_linewidth())
    #     plt.setp(l, linewidth=1)
    #     plt.setp(l, markersize=1)
    # ax.set_xticklabels(["Spontaneous\nActivity (3s)", "On-phasic\nResponse (75ms)",
    #                     "Sustained\nResponse (925ms)", "Off-phasic\nResponse (75ms)",
    #                     "After\nResponse (900ms)"], rotation=45)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylim(-10, 100)
    ax.set_xlabel("")

    ax.legend(bbox_to_anchor=(1.3, 1), fontsize=8, labelspacing=0.5)
    fig.tight_layout()
    fig.savefig("{}_indExpIDFR.png".format(outBase), dpi=300)


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
    spikeTimesDataStacked.rename(columns={0: "Spike Time (ms)",
                                  "level_1": "Interval"}, inplace=True)

    sns.set(style="whitegrid", rc=mplPars)

    fig, ax = plt.subplots(figsize=(7, 5.6))

    sns.violinplot(hue=mdFN['laborState'], y="Spike Time (ms)", x="Interval",
                data=spikeTimesDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=["r", "b"],
                order=spikeTimesCols, split=True, scale="area", inner=None, cut=0, linewidth=0, scale_hue=False)

    sns.pointplot(hue=mdFN['laborState'], y="Spike Time (ms)", x="Interval",
                data=spikeTimesDataStacked, ax=ax, hue_order=['Newly Emerged', 'Forager'], palette=['w', 'w'], join=False,
                  markers="_", dodge=True, ci=None, size=30)

    for colInd, (col, s) in enumerate(spikeTimesData.iteritems()):
        t, pVal = ttest_ind(s["Forager"], s["Newly Emerged"], equal_var=False, nan_policy='omit')
        col = 'k'
        if pVal < 0.05:
            col = 'g'

        ax.text(colInd, -5, "{:1.1e}".format(pVal), fontdict={'color': col}, fontsize=plt.rcParams['xtick.labelsize'],
                horizontalalignment='center', verticalalignment='center')

    xtickLabels = [x[:-5] for x in spikeTimesCols]
    ax.set_xticklabels(xtickLabels, rotation=45)
    ax.set_ylim(-10, 80)
    ax.set_xlabel("")

    # l1, = ax.plot((), (), 'rs', label="Newly Emerged")
    # l2, = ax.plot((), (), 'bs', label="Forager")
    # ax.legend(handles=[l1, l2], loc="upper right")

    ax.legend().set_visible(False)
    
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


    for expInd, (expID, expDF) in enumerate(dataDF.groupby("Experiment ID")):

        print("Doing {}".format(expID))

        expCol = plt.cm.Spectral((expInd + 1) / float(nExpIDs))

        allAx.plot(expDF[spikeFRSpikeTimesFNs["laterFR"]], expDF[spikeFRSpikeTimesFNs["spontFR3"]],
                   color=expCol, ms=5, marker="o", ls="None")

        fig, ax = plt.subplots(figsize=(7, 5.6))
        maxTrials = max(float(s[5:]) + 1 for s in expDF["TrialName"])
        for rowID, rowS in expDF.iterrows():
            trialNo = float(rowS["TrialName"][5:]) + 1
            col = plt.cm.Spectral(1 - trialNo / maxTrials)
            ax.plot(rowS[spikeFRSpikeTimesFNs["laterFR"]],
                    rowS[spikeFRSpikeTimesFNs["spontFR3"]],
                    color=col, marker="o", ms=10)

        ax.set_xlabel("FR during Sustained \n Response (75 to 1000ms)")
        ax.set_ylabel("FR during Spontaneous \n Response (-3000 to 0ms)")
        ax.set_title("VIBGYOR colormap\nViolet=Trial {}, Red=Trial {}".format(1, maxTrials))

        ax.set_xlim(-5, 50)
        ax.set_ylim(-5, 50)

        fig.tight_layout()
        fig.savefig(os.path.join(outDir, "{}.png".format(expID)), dpi=300)
        plt.close(fig.number)

    allAx.set_xlabel("FR during Sustained \n Response (75 to 1000ms)")
    allAx.set_ylabel("FR during Spontaneous \n Response (-3000 to 0ms)")
    maxX = 25
    allAx.plot((0, maxX), (0, maxX), 'k:')
    allAx.set_xlim(-5, 50)
    allAx.set_ylim(-5, allAx.get_ylim()[1] + 10)

    allFig.tight_layout()
    allFig.savefig(os.path.join(outDir, "allExps.png"), dpi=300)

    # allAx.set_ylim(-5, 15)
    # allAx.set_xlim(-5, 20)
    #
    # # allFig.tight_layout()
    # allFig.savefig(os.path.join(outDir, "allExps_zoomed.png"), dpi=300)

    mplPars["legend.frameon"] = True
    mplPars["legend.framealpha"] = 1
    mplPars["legend.fontsize"] = 18
    sns.set(style="whitegrid", rc=mplPars)

    fvneFig, fvneAx = plt.subplots(figsize=(7, 5.6))
    fvneFigMeans, fvneAxMeans = plt.subplots(figsize=(7, 5.6))

    def add_regression_plot(dataDF, x, y, color, ax1):

        xData, yData = dataDF[x], dataDF[y]
        xMin, xMax = xData.min(), xData.max()

        LRRes = linregress(xData, yData)
        ax1.plot(xData, yData, color=color, marker='o', ls='None', ms=5)
        ax1.plot((xMin, xMax), (xMin * LRRes[0] + LRRes[1],
                                                   xMax * LRRes[0] + LRRes[1]), color=color, ls='-', marker='None')
        l1, = ax1.plot((), (), color=color, marker='o', ls='-', ms=5)
        return l1, LRRes

    foragerData = dataDF.loc[lambda x: x[mdFN["laborState"]] == "Forager", :]

    fLegendHandle, fLRRes = add_regression_plot(foragerData, y=spikeFRSpikeTimesFNs["spontFR3"],
                                        x=spikeFRSpikeTimesFNs["laterFR"], color=(0, 0, 1, 1), ax1=fvneAx)

    neData = dataDF.loc[lambda x: x[mdFN["laborState"]] == "Newly Emerged", :]
    neLegendHandle, neLRRes = add_regression_plot(neData, y=spikeFRSpikeTimesFNs["spontFR3"],
                                        x=spikeFRSpikeTimesFNs["laterFR"], color=(1, 0, 0, 1), ax1=fvneAx)

    fvneAx.set_xlabel("FR during Sustained \n Response (75 to 1000ms)")
    fvneAx.set_ylabel("FR during Spontaneous \n Response (-3000 to 0ms)")

    fvneAx.set_xlim(-5, 50)
    fvneAx.set_ylim(-5, fvneAx.get_ylim()[1] + 10)

    fvneAx.legend((fLegendHandle, neLegendHandle),
                  ("Forager,\nslope={:0.3g}".format(fLRRes[0]),
                   "Newly Emerged,\nslope={:0.3g}".format(neLRRes[0])),
                  loc="upper right")
    fvneFig.tight_layout()
    fvneFig.savefig(os.path.join(outDir, "allExpsFVsNE.png"), dpi=300)

    foragerDataMeans = foragerData.groupby(mdFN["expID"]).mean()
    neDataMeans = neData.groupby(mdFN["expID"]).mean()
    fLegendHandle, fLRResMeans = add_regression_plot(foragerDataMeans, y=spikeFRSpikeTimesFNs["spontFR3"],
                                                x=spikeFRSpikeTimesFNs["laterFR"], color=(0, 0, 1, 1), ax1=fvneAxMeans)

    neLegendHandle, neLRResMeans = add_regression_plot(neDataMeans, y=spikeFRSpikeTimesFNs["spontFR3"],
                                                  x=spikeFRSpikeTimesFNs["laterFR"], color=(1, 0, 0, 1), ax1=fvneAxMeans)
    # Ref: https://stats.stackexchange.com/questions/93540/testing-equality-of-coefficients-from-two-different-regressions#99536
    slopeComparisonStat = (fLRResMeans[0] - neLRResMeans[0]) / (np.linalg.norm([fLRResMeans[4], neLRResMeans[4]]))
    nForagerMeans = foragerDataMeans.shape[0]
    nNEMeans = neDataMeans.shape[0]
    foragerMeansVar = (fLRResMeans[4] ** 2) * nForagerMeans
    NEMeansVar = (neLRResMeans[4] ** 2) * nNEMeans
    welchDF = getWelchDF(foragerMeansVar, NEMeansVar, nForagerMeans, nNEMeans, nForagerMeans, nNEMeans)

    pValSlopeDiff = 2 * (1-t.cdf(slopeComparisonStat, welchDF))

    print("Significance of difference of slopes for means: {}".format(pValSlopeDiff))

    fvneAxMeans.set_xlabel("Firing Rate during\nInhibitory Response (spikes/s)")
    fvneAxMeans.set_ylabel("Firing Rate during\nSpontaneous Activity\n(spikes/s)")

    fvneAxMeans.plot((0, 10), (0, 10), 'k--')

    fvneAxMeans.set_xlim(0, 22)
    fvneAxMeans.set_ylim(0, 22)

    fvneAxMeans.legend((neLegendHandle, fLegendHandle),
                       ("Newly\nEmerged,\nslope={:0.3g}".format(neLRResMeans[0]),
                       "Forager,\nslope={:0.3g}".format(fLRResMeans[0])),
                       loc="best")

    # fvneAxMeans.set_aspect("equal")
    fvneFigMeans.tight_layout()
    fvneFigMeans.savefig(os.path.join(outDir, "allExpsFVsNE_means.png"), dpi=300)


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


def plotRelativeInhibition(dataXL, outBase):
    """
    This function calculated relative inhibition as the ratio of firing rates during inhibition and spontaneous activity
    for every response and compares the distribution between foragers and newly emerged adults.
    :param dataXL: string, path of the excel file generated using the function 'saveData' above.
    :param outBase: string, the generated plot will be saved to the file "<outBase>.png"
    :return:
    """

    sns.set(style="whitegrid", rc=mplPars)

    dataDF = pd.read_excel(dataXL)
    dataDFPositiveSpont = dataDF[dataDF[spikeFRSpikeTimesFNs["spontFR3"]] > 0]
    inhFR = dataDFPositiveSpont[spikeFRSpikeTimesFNs["laterFR"]]
    spontFR3 = dataDFPositiveSpont[spikeFRSpikeTimesFNs["spontFR3"]]
    dataDFPositiveSpont["relativeInhibition"] = 1 - (inhFR / spontFR3)
    dataDFPositiveSpont.loc[:, "temp"] = 0


    foragerData = dataDFPositiveSpont.loc[dataDFPositiveSpont[mdFN["laborState"]] == "Forager", "relativeInhibition"]
    neData = dataDFPositiveSpont.loc[dataDFPositiveSpont[mdFN["laborState"]] == "Newly Emerged", "relativeInhibition"]

    t, pVal = ttest_ind(foragerData, neData, equal_var=False)
    print("ForagerMean={}, ForagerStd={}\nneMean={}, neStd={}\npVal={}".format(foragerData.mean(), foragerData.std(),
                                                                               neData.mean(), neData.std(), pVal))

    fig, ax = plt.subplots(figsize=(7, 5.6))

    sns.violinplot(x="relativeInhibition", y="temp", data=dataDFPositiveSpont,
                   hue="Labor State", split=True, scale="area",
                   inner="quartile", palette=["r", "b"], hue_order=["Newly Emerged", "Forager"], ax=ax, orient='h',
                   bw=0.001)
    ax.set_yticklabels([])
    ax.set_ylabel("Relative frequency\n of occurance")
    ax.set_xlabel("Relative Inhibtion")
    ax.legend(ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig("{}.png".format(outBase), dpi=300)



if __name__ == "__main__":

    assert len(sys.argv) == 4, "Improper Usage! Please use as\n" \
                               "python {currFile} save <input excel file> " \
                               "<output data file>  or \n" \
                               "python {currFile} plotFRFvNE <input Data excel file> " \
                               "<output Base> or\n" \
                               "python {currFile} plotFRIndNrns <input Data excel file> " \
                               "<output Base> or\n" \
                               "python {currFile} plotSpikeTimes <input Data excel file> " \
                               "<output Base> or \n" \
                               "python {currFile} plotLaterFRVsSpontFR <input Data excel file> " \
                               "<output Directory> or \n" \
                               "python {currFile} plotFRRelative2SpontFR <input Data excel file>" \
                               "<output Base>or \n" \
                               "python {currFile} plotRelativeInhibition <input Data excel file>" \
                               "<output Base>".format(currFile=sys.argv[0])

    if sys.argv[1] == "save":
        saveData(*sys.argv[2:])
    elif sys.argv[1] == "plotFRFvNE":
        plotFRFVsNE(*sys.argv[2:])
    elif sys.argv[1] == "plotFRIndNrns":
        plotFRIndNrns(*sys.argv[2:])
    elif sys.argv[1] == "plotSpikeTimes":
        plotSpikeTimes(*sys.argv[2:])
    elif sys.argv[1] == "plotLaterFRVsSpontFR":
        plotInhFRVsSpontFR(*sys.argv[2:])
    elif sys.argv[1] == "plotFRRelative2SpontFR":
        plotFRRelative2SpontFR(*sys.argv[2:])
    elif sys.argv[1] == "plotRelativeInhibition":
        plotRelativeInhibition(*sys.argv[2:])
    else:
        print("Unknown usage! Please use as\n"
              "python {currFile} save <input excel file> " 
              "<output data file>  or \n" 
              "python {currFile} plotFRFvNE <input Data excel file> " 
              "<output Base> or\n" 
              "python {currFile} plotFRIndNrns <input Data excel file> " 
              "<output Base> or\n" 
              "python {currFile} plotSpikeTimes <input Data excel file> " 
              "<output Base> or \n" 
              "python {currFile} plotLaterFRVsSpontFR <input Data excel file> " 
              "<output Directory> or \n" 
              "python {currFile} plotFRRelative2SpontFR <input Data excel file>" 
              "<output Base>or \n"
              "python {currFile} plotRelativeInhibition <input Data excel file>"
              "<output Base>".format(currFile=sys.argv[0])
            )











