import os
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from GJEMS.misc.pandasFuncs import dfInterHueFunc
import pandas as pd
from GJEphys.pdColumnNameMapCont import mdFN, fFN, newFFN
from scipy.stats import spearmanr
import numpy as np
from GJEMS.folderDefs import homeFolder
from GJEMS.viz.matplotlibRCParams import mplPars

# plt.ion()


sns.set(rc=mplPars)


allDataFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                           'contStimAllData_expanded.xlsx')
resDir = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                      'DOEParmsFVsNE')

if not os.path.isdir(resDir):
    os.mkdir(resDir)

freqs = [265]

expNames = [
    '130313-4Rh',
    '130322-1LY',
    '130326-2Rh',
    '130408-1LY',
    '130425-1Al',
    '130501-2Rh',
    '130523-3LY',
    # '130605-1LY',
    '130605-2LY',
    '130705-1LY',
    # '140424-1LY',
    # '140701-1Al',
    '140813-3Al',
    '140930-1Al',
    '140917-1Al',
    '141030-1Al',
]

dfFull = pd.read_excel(allDataFile, index_col=(0, 1, 2))
dfFullRI = dfFull.reset_index()
criterion = dfFullRI[mdFN['expID']].isin(expNames) & dfFullRI[mdFN['freq']].isin(freqs)
df = dfFullRI[criterion]

allFFN = dict(fFN, **newFFN)

kruskalTestDF = dfInterHueFunc(data=df, pars=allFFN.values(), hue=mdFN['laborState'],
                               outLabels=['s', 'p'], func=kruskal, kwargsDict={'nan_policy': 'omit'})

# figs = {}

for k, v in allFFN.iteritems():

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
    # figs[k] = fig


    with sns.axes_style('darkgrid'):

        sns.boxplot(x=mdFN['laborState'], y=v, whis=np.inf,
                      data=df, ax=ax, order=['newly emerged', 'forager'])
        sns.stripplot(x=mdFN['laborState'], y=v, data=df, ax=ax, order=['newly emerged', 'forager'],
                      jitter=True, color='k', size=5)
        ax.set_title('p={:.6f}'.format(kruskalTestDF.loc[v, 'p']))

    fig.tight_layout()
    fig.canvas.draw()
    fig.savefig(os.path.join(resDir, '{}FVsNE.png'.format(k)), dpi=300)
    plt.close(fig.number)

kruskalTestDF.to_excel(os.path.join(resDir, 'KruskalWallisTestsFVsNE.xlsx'), sheet_name='Sheet1')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# figs['laterFRVsSpont1'] = fig
laborStates = df[mdFN['laborState']].unique()
cols = sns.color_palette('Set2', laborStates.size)
colsDict = dict(zip(laborStates, cols))
for ls, lsDF in df.groupby(mdFN['laborState']):

    lsCorr, lsP = spearmanr(lsDF[allFFN['laterFR']], lsDF[allFFN['spontFR1']])
    with sns.axes_style('darkgrid'):
        sns.regplot(data=lsDF, x=allFFN['spontFR1'], y=allFFN['laterFR'], color=colsDict[ls], ax=ax,
                    ci=None, truncate=True)
        ax.plot([], [], color=colsDict[ls], marker='o', ls='-', label='{}; p={:.3f}, r={:.3f}'.format(ls, lsP, lsCorr))
        ax.plot([-10, 100], [-10, 100], 'k:')
    ax.legend(loc='best')
    ax.set_xlim(-10, 100)
    ax.set_ylim(-10, 100)

fig.tight_layout()
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'laterFRVsSpontFR1_FVsNE.png'), dpi=300)
plt.close(fig.number)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# figs['initFRVsSpont1'] = fig
laborStates = df[mdFN['laborState']].unique()
cols = sns.color_palette('Set2', laborStates.size)
colsDict = dict(zip(laborStates, cols))
for ls, lsDF in df.groupby(mdFN['laborState']):

    lsCorr, lsP = spearmanr(lsDF[allFFN['initFR']], lsDF[allFFN['spontFR1']])
    with sns.axes_style('darkgrid'):
        sns.regplot(data=lsDF, x=allFFN['spontFR1'], y=allFFN['initFR'], color=colsDict[ls], ax=ax,
                    ci=None, truncate=True)
        ax.plot([], [], color=colsDict[ls], marker='o', ls='-', label='{}; p={:.3f}, r={:.3f}'.format(ls, lsP, lsCorr))
        ax.plot([-10, 100], [-10, 100], 'k:')
    ax.legend(loc='best')
    ax.set_xlim(-10, 100)
    ax.set_ylim(-10, 100)

fig.tight_layout()
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'initFRVsSpontFR1_FVsNE.png'), dpi=300)
plt.close(fig.number)


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# figs['reboundFRVsSpont1'] = fig
laborStates = df[mdFN['laborState']].unique()
cols = sns.color_palette('Set2', laborStates.size)
colsDict = dict(zip(laborStates, cols))
for ls, lsDF in df.groupby(mdFN['laborState']):

    lsCorr, lsP = spearmanr(lsDF[allFFN['reboundFR']], lsDF[allFFN['spontFR1']])
    with sns.axes_style('darkgrid'):
        sns.regplot(data=lsDF, x=allFFN['spontFR1'], y=allFFN['reboundFR'], color=colsDict[ls], ax=ax,
                    ci=None, truncate=True)
        ax.plot([], [], color=colsDict[ls], marker='o', ls='-', label='{}; p={:.3f}, r={:.3f}'.format(ls, lsP, lsCorr))
        ax.plot([-10, 100], [-10, 100], 'k:')
    ax.legend(loc='best')
    ax.set_xlim(-10, 100)
    ax.set_ylim(-10, 100)

fig.tight_layout()
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'reboundFRVsSpontFR1_FVsNE.png'), dpi=300)
plt.close(fig.number)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# figs['reboundFRVsSpont3'] = fig
laborStates = df[mdFN['laborState']].unique()
cols = sns.color_palette('Set2', laborStates.size)
colsDict = dict(zip(laborStates, cols))
for ls, lsDF in df.groupby(mdFN['laborState']):

    lsCorr, lsP = spearmanr(lsDF[allFFN['reboundFR']], lsDF[allFFN['spontFR3']])
    with sns.axes_style('darkgrid'):
        sns.regplot(data=lsDF, x=allFFN['spontFR3'], y=allFFN['reboundFR'], color=colsDict[ls], ax=ax,
                    ci=None, truncate=True)
        ax.plot([], [], color=colsDict[ls], marker='o', ls='-', label='{}; p={:.3f}, r={:.3f}'.format(ls, lsP, lsCorr))
        ax.plot([-10, 100], [-10, 100], 'k:')
    ax.legend(loc='best')
    ax.set_xlim(-10, 100)
    ax.set_ylim(-10, 100)

fig.tight_layout()
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'reboundFRVsSpontFR3_FVsNE.png'), dpi=300)
plt.close(fig.number)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# figs['afterReboundFRVsSpont1'] = fig
laborStates = df[mdFN['laborState']].unique()
cols = sns.color_palette('Set2', laborStates.size)
colsDict = dict(zip(laborStates, cols))
for ls, lsDF in df.groupby(mdFN['laborState']):

    lsCorr, lsP = spearmanr(lsDF[allFFN['afterReboundFR']], lsDF[allFFN['spontFR1']])
    with sns.axes_style('darkgrid'):
        sns.regplot(data=lsDF, x=allFFN['spontFR1'], y=allFFN['afterReboundFR'], color=colsDict[ls], ax=ax,
                    ci=None, truncate=True)
        ax.plot([], [], color=colsDict[ls], marker='o', ls='-', label='{}; p={:.3f}, r={:.3f}'.format(ls, lsP, lsCorr))
        ax.plot([-10, 100], [-10, 100], 'k:')
    ax.legend(loc='best')
    ax.set_xlim(-10, 100)
    ax.set_ylim(-10, 100)

fig.tight_layout()
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'afterReboundFRVsSpontFR1_FVsNE.png'), dpi=300)
plt.close(fig.number)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# figs['reboundFRVsAfterReboundFR'] = fig
laborStates = df[mdFN['laborState']].unique()
cols = sns.color_palette('Set2', laborStates.size)
colsDict = dict(zip(laborStates, cols))
for ls, lsDF in df.groupby(mdFN['laborState']):

    lsCorr, lsP = spearmanr(lsDF[allFFN['afterReboundFR']], lsDF[allFFN['reboundFR']])
    with sns.axes_style('darkgrid'):
        sns.regplot(data=lsDF, x=allFFN['afterReboundFR'], y=allFFN['reboundFR'], color=colsDict[ls], ax=ax,
                    ci=None, truncate=True)
        ax.plot([], [], color=colsDict[ls], marker='o', ls='-', label='{}; p={:.3f}, r={:.3f}'.format(ls, lsP, lsCorr))
        ax.plot([-10, 100], [-10, 100], 'k:')
    ax.legend(loc='best')
    ax.set_xlim(-10, 100)
    ax.set_ylim(-10, 100)

fig.tight_layout()
fig.canvas.draw()
fig.savefig(os.path.join(resDir, 'reboundFRVsAfterReboundFR_FVsNE.png'), dpi=300)
plt.close(fig.number)

