from GJEMS.ephys.simParamExtract import *
from easygui import fileopenbox
import matplotlib.pyplot as plt
import numpy as np

hdf5 = fileopenbox('Indicate the HDF5 file containing Current Clamp Data', filetypes=['*.hdf5'])

spe = SimulationParameterExtracter()

spe.loadCCData(hdf5)
# spe.vizSingleAmpResp()
spe.vizRinAndTaum()
# spe.plotFittedExp()