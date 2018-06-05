import os
from matplotlib import pyplot as plt
import seaborn as sns
from GJEMS.misc.pandasFuncs import dfInterHueFunc
from scipy.stats import kruskal, spearmanr
import pandas as pd
from GJEphys.pdColumnNameMapCont import mdFN, fFN, newFFN
from scipy.stats import spearmanr
import numpy as np
from matplotlib.ticker import MultipleLocator

plt.ion()
mplPars = {'text.usetex': True,
           'axes.labelsize': 'large',
           'axes.titlesize': 42,
           'font.family': 'sans-serif',
           'font.sans-serif': 'computer modern roman',
           'font.size': 42,
           'font.weight': 'black',
           'xtick.labelsize': 36,
           'ytick.labelsize': 36,
           'legend.fontsize': 30,
           }

sns.set_style('darkgrid')
sns.set(rc=mplPars)

homeFolder = os.path.expanduser('~')
allDataFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                           'contStimAllData_expanded.xlsx')
resDir = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'OPExcitationFitting2Exp',
                      'SpikeAdaptation')

if not os.path.isdir(resDir):
    os.mkdir(resDir)

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
freqs = [100, 200, 265, 300, 400]

dfFull = pd.read_excel(allDataFile, index_col=(0, 1, 2))
dfFullRI = dfFull.reset_index()
nrnCrit = dfFullRI[mdFN['expID']].isin(expNames)
crit1307051LY = dfFullRI[mdFN['expID']].isin(['130705-1LY'])
critForager = dfFullRI[mdFN['laborState']].isin(['forager'])
crit265 = dfFullRI[mdFN['freq']].isin(freqs)
df265 = dfFullRI[nrnCrit & crit265]
dfAllFreq = dfFullRI[crit1307051LY]
# dfAllFreq = dfFullRI[critForager]
allFFN = dict(fFN, **newFFN)

latencies = [allFFN['firstSpikeLatency'], allFFN['secondSpikeBISI'],
          allFFN['thirdSpikeBISI'], allFFN['fourthSpikeBISI']]
fig1, axs = plt.subplots(nrows=len(latencies), ncols=1, figsize=(14, 11.2), sharex='all', sharey='all')
cols = sns.color_palette('Set2', n_colors=len(latencies))

for ind, l in enumerate(latencies):
    ax = axs[ind]
    ax.plot([0, 15], [0, 15], color='k', marker='None', ls='--')
    r, p = spearmanr(1000 / dfAllFreq[mdFN['freq']], dfAllFreq[l], nan_policy='omit')
    sns.regplot(x=1000/dfAllFreq[mdFN['freq']], y=dfAllFreq[l],
                color=cols[ind], ax=ax, ci=0, truncate=True,
                label='{}\nr={:1.2f}, p={:1.3f}'.format(l, np.round(r, 2), np.round(p, 3)))
    ax.set_xlim(0, 15)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend(loc='upper left')
axs[-1].set_xlabel('Sine stimulus Period (in ms)')
fig1.tight_layout()
# fig.savefig(os.path.join(resDir, 'latenciesVsFreq.png'), dpi=300)
# plt.close(fig.number)

# mdSub = [mdFN['freq']]
# fnSub = [allFFN['secondSpikeBISI'], allFFN['thirdSpikeBISI'], allFFN['fourthSpikeBISI']]
# dfSub = dfAllFreq.loc[:, mdSub + fnSub]
# dfSub.set_index(mdSub, inplace=True)
# dfSub1 = dfSub.stack()
# dfSub1.name = 'value'
# dfSub2 = dfSub1.reset_index()
# dfSub2.rename(columns={'level_1': 'ISI/Latency'}, inplace=True)
# fig2, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
# sns.boxplot(data=dfSub2, x=mdFN['freq'], y='value', hue='ISI/Latency', ax=ax)
# ax.set_ylim(0, 30)
# ax.set_xlabel('')
# ax.set_ylabel('')
# fig2.tight_layout()
# # fig.savefig(os.path.join(resDir, 'ISIAdaptationVsFreq.png'), dpi=300)
# # plt.close(fig.number)



mdSub = [mdFN['laborState']]
fnSub = [allFFN['firstSpikeLatency'], allFFN['secondSpikeBISI'],
          allFFN['thirdSpikeBISI'], allFFN['fourthSpikeBISI']]
dfSub = df265.loc[:, mdSub + fnSub]
dfSub.set_index(mdSub, inplace=True)
dfSub1 = dfSub.stack()
dfSub1.name = 'value'
dfSub2 = dfSub1.reset_index()
dfSub2.rename(columns={'level_1': 'ISI/Latency'}, inplace=True)

kruskalTestDF = dfInterHueFunc(data=dfSub.reset_index(), pars=fnSub, hue=mdFN['laborState'],
                               outLabels=['s', 'p'], func=kruskal, kwargsDict={'nan_policy': 'omit'})
fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
sns.boxplot(data=dfSub2, x='ISI/Latency', y='value', hue=mdFN['laborState'], ax=ax)
for rowInd, (th, row) in enumerate(kruskalTestDF.iterrows()):
    p = row['p']
    if p < 0.01:
        toPrint = r'p\textless 1\%'
    elif p < 0.05:
        toPrint = r'p\textless 5\%'
    else:
        toPrint = r'p\textgreater 5\%'
    ax.text(rowInd, -5, toPrint, ha='center', va='center')
ax.set_ylim(-10, 40)
ax.set_xlabel('')
ax.set_ylabel('')
ax.legend(loc='upper left')
fig3.tight_layout()

fig4, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 11.2))
tMax = 75  # in ms
xt = np.arange(0, tMax, 0.01)  # in ms
stim = np.sin(2 * np.pi * 265 * 1e-3 * xt - np.pi )

def stimTime2Phase(xt):
    return 2 * np.pi * 265 * 1e-3 * xt - np.pi

def stimPhase2Time(ph):
    return (ph + np.pi) / (2 * np.pi * 265 * 1e-3)

ax.plot(xt, stim, 'r:', label='stimulus')
cols = sns.color_palette(n_colors=2)
dfSub['Second Spike\nLatency (ms)'] = dfSub.apply(lambda x: x[allFFN['firstSpikeLatency']] + x[allFFN['secondSpikeBISI']],
                                                  axis=1)
dfSub['Third Spike\nLatency (ms)'] = dfSub.apply(lambda x: x['Second Spike\nLatency (ms)'] +
                                                           x[allFFN['thirdSpikeBISI']],
                                                  axis=1)
dfSub['Fourth Spike\nLatency (ms)'] = dfSub.apply(lambda x: x['Third Spike\nLatency (ms)'] +
                                                            x[allFFN['fourthSpikeBISI']],
                                                  axis=1)
spLatencies = [allFFN['firstSpikeLatency'], 'Second Spike\nLatency (ms)',
               'Third Spike\nLatency (ms)', 'Fourth Spike\nLatency (ms)']
foragerYs = [1, 0.75, 0.5, 0.25]
neYs = [-0.25, -0.5, -0.75, -1]
dfSub2F = dfSub2.loc[lambda x: x[mdFN['laborState']] == 'forager', :]
dfSub2NE = dfSub2.loc[lambda x: x[mdFN['laborState']] == 'newly emerged', :]

for ind, spl in enumerate(spLatencies):

    splForagerVals = dfSub.loc['forager', spl]
    splForagerVals = splForagerVals[~pd.isnull(splForagerVals)]

    splNEVals = dfSub.loc['newly emerged', spl]
    splNEVals = splNEVals[~pd.isnull(splNEVals)]

    ax.plot(splForagerVals, [foragerYs[ind]] * len(splForagerVals),
            color=cols[0], marker='^', ls='None', ms=10)
    ax.plot(splNEVals, [neYs[ind]] * len(splNEVals),
            color=cols[1], marker='^', ls='None', ms=10)
    ax.boxplot(x=[splForagerVals], vert=False, positions=[foragerYs[ind]],
               boxprops={'color': cols[0]})
    ax.boxplot(x=[splNEVals], vert=False, positions=[neYs[ind]],
               boxprops={'color': cols[1]})

ax.plot(0, 0, color=cols[0], marker='^', ls='None', label='forager', ms=10)
ax.plot(0, 0, color=cols[1], marker='^', ls='None', label='newly emerged', ms=10)
ax.legend(loc='upper right')
ax.set_ylim(-1.25, 2)
ax.set_yticks(foragerYs + neYs)
ax.set_yticklabels(['$1^{st}$', '$2^{nd}$', '$3^{rd}$', '$4^{th}$'] * 2)
ax.set_ylabel('Spike Number')
ax.set_xlabel('Stimulus Phase (radians)')
xticksPhase = np.arange(0, stimTime2Phase(tMax), np.pi)
ax.set_xticks(stimPhase2Time(xticksPhase))
ax.set_xticklabels([r'{:d}$\pi$'.format(int(x / np.pi)) for x in xticksPhase], rotation=90)
# http://matplotlib.org/examples/pylab_examples/major_minor_demo1.html
# minorLocator = MultipleLocator(np.pi / 4)
# ax.xaxis.set_minor_locator(minorLocator)
ax.grid(True, which='both', axis='x')
fig4.tight_layout()
