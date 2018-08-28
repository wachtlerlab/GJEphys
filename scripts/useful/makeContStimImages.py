'''
Description: This script is used to plot time traces of every stimulus and its corresponding response for every
<Experiment ID> specified and saving the figures generated into image files. The files are arranged in a hierarchical
directory structure with neuron category at the first level, <Experiment ID> at the second level and frequency of
stimulus at the third level. Directories are created if they don't yet exit. Existing image files are overwritten.
In addition, an overview figure is generated with all the responses stacked and saved into a file.
<Experiment ID>s can be specified in one of the following two ways:
1. Using the first comment block below to specify a list of neuron categories in the variable "categories". All
<Experiment ID>s of all specified neuron categories are used for generating plots.
2. Using the second or third comment block below to specify a dictionary "expIDsByCat" which contains neuron categories
as keys and list of corresponding <Experiment ID>s as values. See examples in the second and third comment blocks
below.

The work horse of this script is the function GJEphys.expNamePlots.makeContImages. See the documentation of this
function for more info.

Usage: python <path to parent directory>/makeContStimImages.py
'''

import os
from GJEphys.expNamePlots import makeContImages
from GJEphys.folderDefs import homeFolder, NIXPath
import json
import tempfile
from GJEphys.matplotlibRCParams import mplPars

mplPars["lines.linewidth"] = 0.75
mplPars["lines.markersize"] = 4

# **********************************************************************************************************************
from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory
from GJEphys.folderDefs import NIXPath, spike2Path, excelSheet, excel

categories = [
                # 'DL-Int-1',
                # 'DL-Int-2',
                'JO neuron',
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

# expIDsByCat = {'DL-Int-1': ['130523-3LY']}
# # expIDsByCat = {'JO terminal local neuron': ["130415-1LY"]}
# # expIDsByCat = {"JO neuron": ["130517-1Al"]}
# # expIDsByCat = {"JO neuron": ["130318-3LY"]}
# # expIDsByCat = {"JO neuron": ["130704-1LY"]}
# expIDsByCat = {"BilateralN": ["130529-2Al",
#             "131113-1Al"]}
# # expIDsByCat = {"DL-Int-1": ["121217-1"]}
# **********************************************************************************************************************

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


resDir = os.path.join(homeFolder, 'DataAndResults/ephys/ContStimsCatalogue')
if not os.path.exists(resDir):
    os.mkdir(resDir)

type2Color = {'BeforeStimulus': 'r', 'DuringStimulus': 'g', 'AfterStimulus': 'b'}

downSampleFactor = 7

freqs = None

thrash, parFile = tempfile.mkstemp(dir=resDir, suffix='.json')

for cat, expNames in expIDsByCat.iteritems():

    print('Doing {}'.format(cat))
    catResDir = os.path.join(resDir, cat)
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

        makeContImages(parFile)


os.remove(parFile)



