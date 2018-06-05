from GJEphys.NEOFuncs import sliceAnalogSignal
from GJEphys.rawDataProcess import RawDataProcessor
from matplotlib import pyplot as plt
import seaborn as sns
from GJEMS.folderDefs import homeFolder
import os
import quantities as qu
import numpy as np

mplPars = {'text.usetex': True,
           'axes.labelsize': 'large',
           'axes.titlesize': 24,
           'font.family': 'sans-serif',
           'font.sans-serif': 'computer modern roman',
           'font.size': 24,
           'font.weight': 'black',
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'legend.fontsize': 20,
           'legend.frameon': True,
           'legend.framealpha': 0,
           'legend.fancybox': True,
           'backend': 'Qt4Agg'
           }

sns.set(rc=mplPars)
dirPath = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'NIXFiles')

# expName = '130425-1Al'
# starts = [205.670976, 272.0004, 276.160416 ,
#           280.270368, 286.880496, 289.960272,
#           348.714912]
# delay = 0.075

# expName = '130326-2LY'
# starts = [219.965424, 273.291312]
# delay = 0.35 * qu.s

# expName = '130529-2Al'
# starts = [149.595936, 153.737472, 303.312576]
# delay = 0.05 * qu.s

expName = '130205-2LY'
starts = [304.62768, 308.047824, 315.36792,
          318.897936, 322.65912]
delay = 0.095 * qu.s

rdp = RawDataProcessor(expName=expName, dirpath=dirPath, readOnly=True)

starts = np.array(starts) * qu.s
interval = 0.1 * qu.s



plt.ion()


for startInd, start in enumerate(starts):
    fig, ax = plt.subplots(figsize=(14, 11.2))
    resp = sliceAnalogSignal(rdp.voltageSignal, start + delay,
                             start + interval + delay)
    stim = sliceAnalogSignal(rdp.vibrationSignal, start + delay,
                             start + interval + delay)

    tVec = resp.times - resp.t_start
    tVec.unit = qu.ms
    ax.plot(tVec.magnitude, (2 * startInd + 1) * 100 + resp.magnitude, 'b-')
    ax.plot(tVec.magnitude, 2 * startInd * 100 + stim.magnitude, 'r-')

    ax.set_xlabel('time (s)')
    ax.set_title(start + delay)
    plt.show(block=False)
    fig.tight_layout()
    fig.canvas.draw()

