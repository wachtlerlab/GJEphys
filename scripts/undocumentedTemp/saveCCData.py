from GJEphys.simParamExtract import *
from easygui import fileopenbox
import matplotlib.pyplot as plt
import numpy as np


smr = fileopenbox('Indicate the smr file to use', filetypes=['*.smr'])

spe = SimulationParameterExtracter()

#Values for 13113-1Al
# spe.parseSpike2Data(smr, [1670, 1748])
spe.parseSpike2Data(smr)

spe.extractResponses()

spe.saveCCData()
# plt.plot(spe.currentAmps, spe.timeConstants, 'ro', mfc='r', ms=5)
# spe.plotFittedExp()
# plt.figure()
# tVec = spe.voltageSignal.t_start + np.arange(len(spe.voltageSignal)) * spe.voltageSignal.sampling_period
# plt.show(block=False)

