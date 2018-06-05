import numpy as np
import os
import json
from matplotlib import pyplot as plt
plt.ion()

def accept_test1(x_new):
    [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = x_new

    return all(x > 0 for x in [taur1, taud1, taur2, taud2, A1, A2]) and bool(abs(offset) < 3) \
            and bool(-50 <= onset1 <= 50) and bool(0 <= onset2 <= 100) \
            and all(x < 100 for x in [taur1, taud1, taur2]) and bool(taud2 < 1000)


expNames = [
            '130313-4Rh',
            # '130322-1LY',
            # '130326-2Rh',
            # '130408-1LY',
            # '130425-1Al',
            # '130501-2Rh',
            # '130523-3LY',
            # '130605-1LY',
            # '130605-2LY',
            # '130705-1LY',
            # '140424-1LY',
            # '140701-1Al',
            # '140813-3Al',
            # '140930-1Al',
            # '140917-1Al',
            # '141030-1Al',
            ]

homeFolder = os.path.expanduser('~')
# resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Backups/OPExcitationFitting_2exp_1/')
resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/OPExcitationFitting2Exp')

fig0, ax0 = plt.subplots(figsize=(10, 8))
fig1, ax1 = plt.subplots(figsize=(10, 8))
fig2, ax2 = plt.subplots(figsize=(10, 8))
fig3, ax3 = plt.subplots(figsize=(10, 8))

figs = [fig0, fig1, fig2, fig3]
axs = [ax0, ax1, ax2, ax3]


for expInd, expName in enumerate(expNames):

    print('Doing: ' + expName)

    nrnResDir = os.path.join(resDir, expName)
    trialFiles = sorted([x for x in os.listdir(nrnResDir) if x.endswith(".json")])

    for trialInd, trialFile in enumerate(trialFiles):

        outFile = os.path.join(nrnResDir, trialFile)
        with open(outFile, 'r') as fle:
            pars = json.load(fle)
            p0s = pars['p0s']
            xs = pars['xs']
            fs = pars['fs']

        xs_accepted = []
        fs_accepted = []
        acceptance = []

        for xInd, x in enumerate(xs):

            f = fs[xInd]

            [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset] = x

            if taud1 < taur1:
                [A1, taud1, taur1] = [-A1, taur1, taud1]

            if taud2 < taur2:
                [A2, taud2, taur2] = [-A2, taur2, taud2]

            if A1 < 0 and A2 < 0:
                x = [-A2, -A1, t2, t1, taur2, taud2, taur1, taud1, offset]
            else:
                x = [A1, A2, t1, t2, taur1, taud1, taur2, taud2, offset]

            xs[xInd] = x

            acc = accept_test1(x)
            acceptance.append(acc)

            if acc:
                xs_accepted.append(x)
                fs_accepted.append(f)

        maxF = max(fs_accepted)
        minF = min(fs_accepted)
        # normFs = (np.array(fs_accepted) - minF) / (maxF - minF)

        [ax.clear() for ax in axs]

        sc0 = ax0.scatter([x[0] for x in xs_accepted], [x[1] for x in xs_accepted],
                         c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
        sc1 = ax1.scatter([x[2] for x in xs_accepted], [x[3] for x in xs_accepted],
                            c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
        sc2 = ax2.scatter([x[4] for x in xs_accepted], [x[6] for x in xs_accepted],
                            c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
        sc3 = ax3.scatter([x[5] for x in xs_accepted], [x[7] for x in xs_accepted],
                          c=fs_accepted, cmap=plt.cm.jet, marker='o', vmin=minF, vmax=maxF, edgecolors='face')
        print(sum(acceptance) / float(len(acceptance)))

        ax0.set_xlabel('A1')
        ax0.set_ylabel('A2')
        ax0.axis('square')
        # ax0.set_xlim(-10, 100)
        # ax0.set_ylim(-10, 100)

        ax1.set_xlabel('onset1')
        ax1.set_ylabel('onset2')
        ax1.axis('square')
        ax1.set_xlim(-50, 100)
        ax1.set_ylim(-50, 100)

        ax2.set_xlabel('taur1')
        ax2.set_ylabel('taur2')
        ax2.axis('square')
        # xlim = [0, 1000]
        # ax2.plot(xlim, [x - 50 for x in xlim], '-m')
        # ax2.plot(xlim, [x + 50 for x in xlim], '-m')
        # ax2.set_xlim(-10, 1000)
        # ax2.set_ylim(-10, 1000)

        ax3.set_xlabel('taud1')
        ax3.set_ylabel('taud2')
        ax3.axis('square')


        # fig0.colorbar(sc0, ax=ax0)
        # fig1.colorbar(sc1, ax=ax1)
        # fig2.colorbar(sc2, ax=ax2)
        # fig3.colorbar(sc3, ax=ax3)

        for fig in [fig0, fig1, fig2, fig3]:
            fig.tight_layout()
            fig.canvas.draw()
            fig.axes[0].set_title(expName)

        import ipdb
        ipdb.set_trace()


