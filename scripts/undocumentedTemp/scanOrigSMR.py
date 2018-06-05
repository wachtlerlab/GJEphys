from GJEphys.rawDataImport import RawDataViewer
import matplotlib.pyplot as plt
plt.ion()
import quantities as qu
import numpy as np
import sys

mplPars = {'text.usetex': True,
           'axes.labelsize': 'large',
           'axes.titlesize': 42,
           'font.family': 'sans-serif',
           'font.sans-serif': 'computer modern roman',
           'font.size': 42,
           'font.weight': 'black',
           'xtick.labelsize': 36,
           'ytick.labelsize': 36,
           'legend.fontsize': 20,
           }

plt.rcParams.update(mplPars)

assert len(sys.argv) == 2, 'Improper Usage! Please use as follows:\npython scanOrigSMR.py <smr file>'
smrFile = sys.argv[1]
see = RawDataViewer(smrFile, forceUnits=True)
tStart = see.vibrationSignal.t_start
tStart.units = qu.s
tStop = see.vibrationSignal.t_stop
tStop.units = qu.s

ts = np.arange(tStart.magnitude, tStop.magnitude, 20) * qu.s

fig, ax = plt.subplots(figsize=(14, 11.2))


for ind in range(len(ts) - 1):

    see.plotVibEpoch(ax, [ts[ind], ts[ind + 1]])
    fig.tight_layout()
    raw_input('Press any key to proceed..')
    ax.clear()


