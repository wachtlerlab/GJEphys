
'''
Usage: This script is used to process raw Data in NIX for after importing using 'rawDataImporter'.
The files to be processed could be specified in two exclusive ways:
1. using the first comment block below. All files of experiments of all interneuron categories
will be processed
2. using the second comment block below. All files of experiments with IDs specifed in the list 'expNames' will be
processed.
3. using the third comment block below. All files of experiments of the interneuron categories specified in the list
'categories' will be processed

The crux of this script is a set of functions called in a specific order with the try except block below.
Look at the documentation of these functions for more information.

Note: This script generates an excel file at the path in the variable GJEphys.folderDefs.processedResultsFile. This
excel file contains information whether processing each specified <Experiment ID> succeeded or failed. This excel
file is used in further scripts to determine whether a file has been successfully processed. Further, this also
generates a lot of log output onto stdout, especially about isolated pulses. It might be useful to redirect
this logged information to a files as shown in "Usage" below.

Usage:
"python rawDataProcessor.py > <path to log file>"
'''

from GJEphys.rawDataProcess import RawDataProcessor
from GJEphys.folderDefs import NIXPath, processedResultsFile
import quantities as qu
import os
import pandas as pd
from traceback2 import print_exception
import sys

# **********************************************************************************************************************
# from GJEMS.GJEphys.KKHAXLParsing import parseMetaDataFile
# from GJEMS.folderDefs import excel, excelSheet, spike2Path
# mdDF = parseMetaDataFile(excel, excelSheet, spike2Path)
# expNames = map(str, mdDF.index)
# **********************************************************************************************************************
# expNames = ["130605-1LY"]
# **********************************************************************************************************************
from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory
from GJEphys.folderDefs import excel, excelSheet, spike2Path

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

expIDsByCat = getExpIDsByCategory(excel, excelSheet, categories, spike2Path)
expNames = []
for cat, catExpNames in expIDsByCat.iteritems():
    expNames += catExpNames
# **********************************************************************************************************************


resultsDF = pd.DataFrame()
for expName in expNames:

    nixFile = os.path.join(NIXPath, '{}.h5'.format(expName))
    if not os.path.exists(nixFile):
        print('File not found {}. Skipping it.'.format(nixFile))
        continue
    print('Processing ' + expName)
    try:
        see = RawDataProcessor(expName=expName, dirpath=NIXPath, askShouldReprocess=False)
        see.detectSpikes()
        see.downSampleVibSignal()
        see.getStimulusEpochs(5 * qu.ms)
        see.extractResponse()
        see.sortSegments()
        see.write2nixFile()
        see.close()
        resultsDF = resultsDF.append({"Experiment ID": expName, "Result": True}, ignore_index=True)
    except Exception as exep:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print_exception(type(exep), str(exep), exc_traceback, file=sys.stdout)
        resultsDF = resultsDF.append({"Experiment ID": expName, "Result": False}, ignore_index=True)

resultsDF.set_index("Experiment ID", inplace=True)
resultsDF.to_excel(processedResultsFile)

