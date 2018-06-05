'''
Description: This script is used to generating plots of firing rate vs freq for all <Experiment ID>s specified. One
figure is generated per category containing such plots for all <Experiment IDs>s of that category. In addition, a figure
is generated containing the distributions of firing rates at different frequencies, pooled across <Experiment ID>s and
saved into an image file. These files are arranged into folders, one folder per category.

<Experiment ID>s can be specified in one of the following two ways:
1. Using the first comment block below to specify a list of neuron categories in the variable "categories". All
<Experiment ID>s of all specified neuron categories are used for generating plots.
2. Using the second or third comment block below to specify a dictionary "expIDsByCat" which contains neuron categories
as keys and list of corresponding <Experiment ID>s as values. See examples in the second and third comment blocks
below.

The work horse of this script is the function GJEphys.expNamePlots.makeSpikeTotalVsFreqImages. Look at the documentation
of this function for more info.

Usage:
python <path to parent directory>/makeSpikeTotalVsFreqImages.py
'''

import os
from GJEphys.expNamePlots import makeSpikeTotalVsFreqImages
from GJEphys.folderDefs import homeFolder
import json
import tempfile
from GJEMS.viz.matplotlibRCParams import mplPars


# **********************************************************************************************************************
from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory
from GJEphys.folderDefs import NIXPath, spike2Path, excelSheet, excel

categories = [
                'DL-Int-1',
                'DL-Int-2',
                'JO neuron',
                'MB neurons',
                'BilateralN',
                'DescendingN',
                'AscendingN',
                'Bilateral Descending N',
                'JO terminal local neuron'
            ]

mdDF = parseMetaDataFile(excel, excelSheet, spike2Path)
expNames = map(str, mdDF.index)

expIDsByCat = getExpIDsByCategory(excel, excelSheet, categories, spike2Path)
# **********************************************************************************************************************
# NIXPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
# expIDsByCat = {'DL-Int-1': ['130313-4Rh']}
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


resDir = os.path.join(homeFolder, 'DataAndResults/GJEphys/Results/SpikeRateVsFreq')
if not os.path.exists(resDir):
    os.mkdir(resDir)

freqs = None

thrash, parFile = tempfile.mkstemp(dir=resDir, suffix='.json')

for cat, expNames in expIDsByCat.iteritems():

    expNamesFiltered = []

    for expName in expNames:
        smrFile = os.path.join(NIXPath, "{}.h5".format(expName))
        if os.path.exists(smrFile):
            expNamesFiltered.append(expName)
        else:
            print('File not found {}. Skipping it.'.format(smrFile))

    print('Doing {}: {}'.format(cat, expNamesFiltered))
    catResDir = os.path.join(resDir, cat)
    if not os.path.exists(catResDir):
        os.mkdir(catResDir)

    pars = {"NIXPath": NIXPath,
            "expNamesDict": {cat: expNamesFiltered},
            "freqs": freqs,
            "catResDir": catResDir,
            "mplPars": mplPars}

    with open(parFile, 'w') as fle:

        json.dump(pars, fle)

    makeSpikeTotalVsFreqImages(parFile)

os.remove(parFile)



