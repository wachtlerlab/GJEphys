import os
from GJEphys.pdColumnNameMapCont import mdFN, newFFN, fFN
from GJEphys.pandasFuncs import dfCorr
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr


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
    # '130322-1LY',
    '130326-2Rh',
    '130408-1LY',
    '130425-1Al',
    # '130501-2Rh',
    # '130523-3LY',
    # '130605-1LY',
    # '130605-2LY',
    '130705-1LY',
    # '140424-1LY',
    # '140701-1Al',
    '140813-3Al',
    '140930-1Al',
    '140917-1Al',
    '141030-1Al',
]

freqs = [100, 200, 265, 300, 365]

homeFolder = os.path.expanduser('~')
allDataFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                           'contStimAllData_expanded.xlsx')
resDir = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                      'DOEParmsVsSpontFR')

if not os.path.isdir(resDir):
    os.mkdir(resDir)

dfFull = pd.read_excel(allDataFile, index_col=(0, 1, 2))
dfFullRI = dfFull.reset_index()
criterion = dfFullRI[mdFN['expID']].isin(expNames) & dfFullRI[mdFN['freq']].isin(freqs)
df = dfFullRI[criterion]

allFFN = dict(fFN, **newFFN)

figsFR3 = {}
fr3Dir = os.path.join(resDir, 'SpontFR3')
if not os.path.isdir(fr3Dir):
    os.mkdir(fr3Dir)

for k, v in allFFN.iteritems():
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
    figsFR3[k] = fig

    with sns.axes_style('darkgrid'):
        cols = sns.color_palette('Set2', len(expNames))
        colsDict = {expName: col for expName, col in zip(expNames, cols)}
        for expName, groupDF in df.groupby(mdFN['expID']):
            sns.regplot(data=groupDF, x=allFFN['spontFR3'], y=allFFN[k],
                        color=colsDict[expName], ax=ax, ci=0, order=1, truncate=True, scatter=True)
            ax.plot([], [], color=colsDict[expName], marker='o', ls='None', label=expName)
        ax.legend(loc='best')


    fig.tight_layout()
    fig.canvas.draw()
    fig.savefig(os.path.join(fr3Dir, '{}VsSpontFR3.jpg'.format(k)), dpi=300)
    plt.close(fig.number)

freqCorrs = dfCorr(data=df, x=newFFN['spontFR3'], ys=allFFN.values(), func=spearmanr,
                   outLabels=['r', 'p'], hue=mdFN['expID'])
freqCorrs.to_excel(os.path.join(fr3Dir, 'ParsCorrWithSpontFR3.xlsx'), sheet_name='Sheet1')

figsFR1 = {}
fr1Dir = os.path.join(resDir, 'SpontFR1')
if not os.path.isdir(fr1Dir):
    os.mkdir(fr1Dir)


for k, v in allFFN.iteritems():
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
    figsFR1[k] = fig

    with sns.axes_style('darkgrid'):

        cols = sns.color_palette('Set2', len(expNames))
        colsDict = {expName: col for expName, col in zip(expNames, cols)}
        for expName, groupDF in df.groupby(mdFN['expID']):
            sns.regplot(data=groupDF, x=allFFN['spontFR1'], y=allFFN[k],
                        color=colsDict[expName], ax=ax, ci=0, order=1, truncate=True, scatter=True)
            ax.plot([], [], color=colsDict[expName], marker='o', ls='None', label=expName)
        ax.legend(loc='best')


    fig.tight_layout()
    fig.canvas.draw()
    fig.savefig(os.path.join(fr1Dir, '{}VsSpontFR1.jpg'.format(k)), dpi=300)
    plt.close(fig.number)

freqCorrs = dfCorr(data=df, x=newFFN['spontFR1'], ys=allFFN.values(), func=spearmanr,
                   outLabels=['r', 'p'], hue=mdFN['expID'])
freqCorrs.to_excel(os.path.join(fr1Dir, 'ParsCorrWithSpontFR1.xlsx'), sheet_name='Sheet1')


