from GJEphys.rawDataProcess import *
import matplotlib.pyplot as plt
import quantities as qu
from GJEphys.spikeDetection import SpikeDetector
import json


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

            '130318-3LY',
            '130425-3Al',
            '130517-1Al',
            '130704-1LY',

]

scales = {}
outFile = '/home/ajay/DataAndResults/ephys/scales.json'

for expName in expNames:


    print('Measuring Scale of ' + expName)
    # ipdb.set_trace()
    see = RawDataProcessor(expName=expName, dirpath='/home/ajay/DataAndResults/GJEphys/NIXFiles/', readOnly=True)

    sd = SpikeDetector(see.voltageSignal)
    sd.filterButterworth()

    sd.thresholdAndDetect(0.5, 5)
    scales[expName] = float(sd.spikeHeights.mean())

with open(outFile, 'w') as fle:
    json.dump(scales, fle)