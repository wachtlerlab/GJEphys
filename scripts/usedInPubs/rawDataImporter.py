'''
Description: This script is used to import raw Data in SMR format into NIX. The files to be imported could be specified
in two exclusive ways:
1. using the first comment block below. All files of experiments of the interneuron categories specified in the list
'categories' will be imported
2. using the second comment block below. All files of experiments with IDs specifed in the list 'expNames' will be
imported.
See the documentation of 'importAll' function of 'GJEphys.rawDataImport' for more about the import process.

Usage: python <path to parent directory>/rawDataImporter.py
'''

from GJEphys.rawDataImport import importAll
from GJEphys.folderDefs import NIXPath, importResultsFile
from GJEphys.folderDefs import excel, excelSheet, spike2Path

# **********************************************************************************************************************

# from GJEphys.KKHAXLParsing import getExpIDsByCategory
#
# categories = [
#                 # 'DL-Int-1',
#                 # 'DL-Int-2',
#                 'JO neuron',
#                 # 'MB neurons',
#                 # 'BilateralN',
#                 # 'DescendingN',
#                 # 'AscendingN',
#                 # 'Bilateral Descending N',
#                 # 'JO terminal local neuron'
#             ]
#
# expIDsByCat = getExpIDsByCategory(excel, excelSheet, categories, spike2Path)
# expNames = []
# for cat, catExpNames in expIDsByCat.iteritems():
#     expNames += catExpNames

# **********************************************************************************************************************

expNames = [
            "130214-1Rh",
            # "130430-1Al",
            # "130607-2Al",
            # "140917-1Al",
            # "130415-2Rh",
            # "130605-1LY",
            # "130705-1LY"

            # "130517-1Al",
            # "130318-3LY",
            # "130704-1LY",
            # "130425-3Al"
            ]

# **********************************************************************************************************************

importResultsDF = importAll(smrPath=spike2Path,
          nixPath=NIXPath,
          excelFile=excel,
          excelSheet=excelSheet,
          expNames=expNames,
          askShouldReplace=False,
          forceUnits=True)

# **********************************************************************************************************************



importResultsDF.to_excel(importResultsFile)


