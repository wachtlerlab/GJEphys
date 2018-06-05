import os
from GJEphys.pdColumnNameMapCont import mdFN, fFN, newFFN
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import pandas as pd
from GJEMS.misc.pandasFuncs import dfCorr

# plt.ion()

mplPars = {'text.usetex': True,
           'axes.labelsize': 'large',
           'axes.titlesize': 42,
           'font.family': 'sans-serif',
           'font.sans-serif': 'computer modern roman',
           'font.size': 42,
           'font.weight': 'black',
           'xtick.labelsize': 36,
           'ytick.labelsize': 36,
           'legend.fontsize': 36,
           }

sns.set(rc=mplPars)

expNames = [
            '130313-4Rh',
            '130425-1Al',
            '130705-1LY',
            ]

freqs = [100, 200, 265, 300, 400]


homeFolder = os.path.expanduser('~')
allDataFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                           'contStimAllData_expanded.xlsx')
resDir = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                           'DOEParmsVsFreq')

if not os.path.isdir(resDir):
    os.mkdir(resDir)

dfFull = pd.read_excel(allDataFile, index_col=(0, 1, 2))
dfFullRI = dfFull.reset_index()
criterion = dfFullRI[mdFN['expID']].isin(expNames) & dfFullRI[mdFN['freq']].isin(freqs)
df = dfFullRI[criterion]



allFFN = dict(fFN, **newFFN)

figs = {}

for k, v in allFFN.iteritems():

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
    figs[k] = fig


    with sns.axes_style('darkgrid'):

        cols = sns.color_palette('Set2', len(expNames))
        colsDict = {expName: col for expName, col in zip(expNames, cols)}
        for expName, groupDF in df.groupby(mdFN['expID']):
            groupDF1 = groupDF.loc[~pd.isnull(groupDF[allFFN[k]]), :]
            if groupDF1.size:
                sns.regplot(data=groupDF1, x=mdFN['freq'], y=allFFN[k],
                            color=colsDict[expName], ax=ax, ci=0, truncate=True)
                ax.plot([], [], color=colsDict[expName], marker='o', ls='None', label=expName)
        ax.legend(loc='best')
        ax.set_xlim(75, 425)
        # ax.legend(loc=(-0.15, 1), ncol=2)

    fig.tight_layout()
    # fig.tight_layout(rect=(0, 0, 1, 0.85))
    fig.canvas.draw()
    fig.savefig(os.path.join(resDir, '{}VsFreq.png'.format(k)), dpi=300)
    plt.close(fig.number)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
figs['inhiMaxNormedVsExciMaxNormed'] = fig

with sns.axes_style('darkgrid'):

    cols = sns.color_palette('Set2', len(expNames))
    colsDict = {expName: col for expName, col in zip(expNames, cols)}
    for expName, groupDF in df.groupby(mdFN['expID']):
        sns.regplot(data=groupDF, x=allFFN['exciMaxNormed'], y=allFFN['inhiMaxNormed'],
                   color=colsDict[expName], ax=ax, ci=0, truncate=True)
        ax.plot([], [], color=colsDict[expName], marker='o', ls='None', label=expName)
    ax.legend(loc='best')

fig.tight_layout()
# fig.tight_layout(rect=(0, 0, 1, 0.85))
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'inhiMaxNormedVsExciMaxNormed.png'), dpi=300)

freqCorrs = dfCorr(data=df, x=mdFN['freq'], ys=allFFN.values(), func=spearmanr, outLabels=['r', 'p'], hue=mdFN['expID'])
exciInhiCorrs = dfCorr(df, x=allFFN['exciMaxNormed'], ys=[allFFN['inhiMaxNormed']], func=spearmanr,
                       outLabels=['r', 'p'], hue=mdFN['expID'])

freqCorrs.to_excel(os.path.join(resDir, 'ParsCorrWithFreq.xlsx'), sheet_name='Sheet1')
exciInhiCorrs.to_excel(os.path.join(resDir, 'inhiMaxNormedVsExciMaxNormed.xlsx'), sheet_name='Sheet1')



# for k, v in parameterCollection.iteritems():
#     fig, ax = plt.subplots(figsize=(14, 11.2))
#
#     figs[k] = fig
#
#     ax.boxplot(v.values(), positions=v.keys(), bootstrap=True, sym=None, widths=[20] * len(v),
#                patch_artist=patches.Patch(), boxprops=dict(color='b'),
#                medianprops=dict(color='b'), whiskerprops=dict(color='b'))
#     ax.plot(allFreqs, map(np.mean, [v[x] for x in allFreqs]), 'r-')
#
#     allFreqVals = []
#     allParamVals = []
#     for freq, freqParamValues in v.iteritems():
#
#         ax.plot([freq] * len(freqParamValues), freqParamValues, 'bo')
#         allFreqVals += [freq] * len(freqParamValues)
#         allParamVals += freqParamValues
#
#     parSpearmanR, parPVal = spearmanr(allFreqVals, allParamVals)
#
#     ax.set_xlabel('Frequency (Hz)')
#     ax.set_ylabel(k + ' (' + units[k] + ')')
#     ax.set_xlim(min(allFreqs) - 20, max(allFreqs) + 20)
#     ax.set_title('spearmanRVal=' + str(parSpearmanR) + '; spearmanPVal=' + str(parPVal))
#
#     fig.tight_layout()
#     fig.canvas.draw()
#
#
# fig, ax = plt.subplots(figsize=(14, 11.2))
# allExciMaxNormedVals = reduce(operator.add, [parameterCollection['ExciMaxNormed'][freq] for freq in allFreqs])
# allInhiMaxNormedVals = reduce(operator.add, [parameterCollection['InhiMaxNormed'][freq] for freq in allFreqs])
# ax.plot(allExciMaxNormedVals, allInhiMaxNormedVals, 'bo')
# exciInhiR, exciInhiP = spearmanr(allExciMaxNormedVals, allInhiMaxNormedVals)
# ax.set_title('spearmanRVal=' + str(exciInhiR) + '; spearmanPVal=' + str(exciInhiP))
# ax.set_xlabel('ExciMaxNormed')
# ax.set_ylabel('InhiMaxNormed')

