'''
Description: This script is used to make figures containing time traces of pulse train stimuli and their corresponding
responses and save them into image files. The files are arranged in a hierarchical directory structure with neuron \
category at the first level, <Experiment ID> at the second level and the tuple
(Pulse Duration, Pulse interval, stimulus frequency) on the third level. Directories are created if they don't yet exit.
Existing image files are overwritten.
In addition, an overview figure is generated with all the responses stacked and saved into a file.
<Experiment ID>s can be specified in one of the following two ways:
1. Using the first comment block below to specify a list of neuron categories in the variable "categories". All
<Experiment ID>s of all specified neuron categories are used for generating plots.
2. Using the second or third comment block below to specify a dictionary "expIDsByCat" which contains neuron categories
as keys and list of corresponding <Experiment ID>s as values. See examples in the second and third comment blocks
below.

The work horse of this script is the function GJEphys.expNamePlots.makePulseImages. See the documentation of this
function for more info.

Usage: python <path to parent directory>/makePulseStimImages.py

'''
import os
from matplotlib import pyplot as plt
from GJEphys.expNamePlots import makePulseImages
from GJEphys.matplotlibRCParams import mplPars
import json
import tempfile


plt.rcParams.update(mplPars)
mplPars["lines.linewidth"] = 0.75
mplPars["lines.markersize"] = 4

# **********************************************************************************************************************
from GJEphys.folderDefs import homeFolder, excel, excelSheet, NIXPath, spike2Path
from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory

categories = [
                'DL-Int-1',
                # 'DL-Int-2',
                # 'JO neuron',
                # 'MB neurons',
                # 'BilateralN',
                # 'DescendingN',
                # 'AscendingN',
                # 'Bilateral Descending N',
                # 'JO terminal local neuron'
            ]

mdDF = parseMetaDataFile(excel, excelSheet, spike2Path)
expNames = map(str, mdDF.index)

expIDsByCat = getExpIDsByCategory(excel, excelSheet, categories, spike2Path)
# **********************************************************************************************************************
# # NIXPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
# # # expIDsByCat = {'DL-Int-1': ['130313-4Rh']}
# # expIDsByCat = {"BilateralN": ["130529-2Al",
# #             "131113-1Al"]}
# expIDsByCat = {'DL-Int-1': ['130605-1LY']}
# **********************************************************************************************************************
# NIXPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
# expIDsByCat = {
#                 'AscendingN':               ['130205-1LY'],
#                 'Bilateral Descending N':   ['130313-3LY'],
#                 'BilateralN':               ['130514-1LY'],
#                 'DL-Int-1':                 ['140424-1LY'],
#                 'DL-Int-2':                 ['130205-2LY'],
#                 'DescendingN':              ['130320-1Rh'],
#                 'JO neuron':                ['130517-1Al'],
#                 'JO terminal local neuron': ['130206-2Rh'],
#                 'MB neurons':               ['130418-2Al']
#                 }

# **********************************************************************************************************************

outDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/PulseStimsCatalogue')
if not os.path.exists(outDir):
    os.mkdir(outDir)

type2Color = {'BeforeStimulus': 'r', 'AfterStimulus': 'b', 'OnPulses': 'g', 'OffPulses': 'y'}


downSampleFactor = 7
thrash, parFile = tempfile.mkstemp(dir=outDir, suffix='.json')

freqs = None

for cat, expNames in expIDsByCat.iteritems():

    print('Doing {}'.format(cat))
    catResDir = os.path.join(outDir, cat)
    if not os.path.exists(catResDir):
        os.mkdir(catResDir)

    for expName in expNames:

        smrFile = os.path.join(NIXPath, '{}.h5'.format(expName))
        if not os.path.exists(smrFile):
            print('File not found {}. Skipping it.'.format(smrFile))
            continue
        print('Doing {}'.format(expName))

        pars = {'NIXPath': NIXPath,
                'expName': expName,
                'freqs': freqs,
                'catResDir': catResDir,
                'downSampleFactor': downSampleFactor,
                'type2Color': type2Color,
                'mplPars': mplPars}

        with open(parFile, 'w') as fle:

            json.dump(pars, fle)

        makePulseImages(parFile)

os.remove(parFile)
















