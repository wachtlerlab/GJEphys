import numpy as np
import os
import json
from matplotlib import pyplot as plt
plt.ion()

def accept_test1(x_new):
    [A1, onset1, tau1, A2, onset2, tau2, offset] = x_new

    return all(x > 0 for x in [tau1, tau2, A1, A2]) and bool(abs(offset) < 3)


expNames = [
            '130313-4Rh',
            '130322-1LY',
            '130326-2Rh',
            '130408-1LY',
            '130425-1Al',
            '130501-2Rh',
            '130523-3LY',
            '130605-1LY',
            '130605-2LY',
            '130705-1LY',
            '140424-1LY',
            '140701-1Al',
            '140813-3Al',
            '140930-1Al',
            '140917-1Al',
            '141030-1Al',
            ]

homeFolder = os.path.expanduser('~')
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Backups/OPExcitationFitting_alpha_1/')




for expInd, expName in enumerate(expNames):

    print('Doing: ' + expName)
    outFile = os.path.join(resDir, expName + '.json')
    with open(outFile, 'r') as fle:
        pars = json.load(fle)
        p0s = pars['p0s']
        xs = pars['xs']
        fs = pars['fs']

    xs_accepted = []
    fs_accepted = []
    acceptance = []

    for x, f in zip(xs, fs):

        [A1, onset1, tau1, A2, onset2, tau2, offset] = x

        if A1 < 0 and A2 < 0:
            x = [-A2, onset2, tau2, -A1, onset1, tau1, offset]

        acc = accept_test1(x)
        acceptance.append(acc)

        if acc:
            xs_accepted.append(x)
            fs_accepted.append(f)

    maxF = max(fs_accepted)
    minF = min(fs_accepted)
    # normFs = (np.array(fs_accepted) - minF) / (maxF - minF)

    fig0, ax0 = plt.subplots(figsize=(10, 8))
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    sc0 = ax0.scatter([x[0] for x in xs_accepted], [x[3] for x in xs_accepted],
                     c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
    sc1 = ax1.scatter([x[1] for x in xs_accepted], [x[4] for x in xs_accepted],
                        c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
    sc2 = ax2.scatter([x[2] for x in xs_accepted], [x[5] for x in xs_accepted],
                        c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
    print(sum(acceptance) / float(len(acceptance)))

    ax0.set_xlabel('A1')
    ax0.set_ylabel('A2')
    ax0.axis('square')
    ax0.set_xlim(-10, 100)
    ax0.set_ylim(-10, 100)

    ax1.set_xlabel('onset1')
    ax1.set_ylabel('onset2')
    ax1.axis('square')
    # ax1.set_xlim(-50, 100)
    # ax1.set_ylim(-50, 100)

    ax2.set_xlabel('tau1')
    ax2.set_ylabel('tau2')
    ax2.axis('square')
    xlim = [0, 1000]
    ax2.plot(xlim, [x - 50 for x in xlim], '-m')
    ax2.plot(xlim, [x + 50 for x in xlim], '-m')
    # ax2.set_xlim(-10, 1000)
    # ax2.set_ylim(-10, 1000)


    fig0.colorbar(sc0, ax=ax0)
    fig1.colorbar(sc1, ax=ax1)
    fig2.colorbar(sc2, ax=ax2)

    for fig in [fig0, fig1, fig2]:
        fig.tight_layout()
        fig.canvas.draw()
        fig.axes[0].set_title(expName)

    import ipdb
    ipdb.set_trace()

    plt.close('all')
