'''
Description: The script generates a table summarizing the number of trials and number of experiments for every category
of neuron and for every frequency of stimulus applied. This script uses NIX files in the directory "NIXPath"
in GJEphy.folderDefs. The original smr files should be imported using "rawDataImporter.py" and the resulting files
processed using "rawDataProcessor.py" before running this script to generate a table of summary statistics. This
script checks whether the NIX file of an input <Experiment ID> has been successfully processed by looking at the
excel file at the path indicated by GJEphys.folderDefs.processedResultsFile.

Usage: python <path to parent directory>/contStimStats.py
'''
from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory
from GJEphys.rawDataAnalyse import RawDataAnalyser
import pylatex
import pandas as pd
from GJEphys.folderDefs import NIXPath, excel, excelSheet, spike2Path, contStimsResultsFile, \
    processedResultsFile, contStimDetailedStatsFile

processedResultsDF = pd.read_excel(processedResultsFile, index_col=0, header=0)

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
columns = ['Category', 'Experiment ID', 'Frequency (Hz)', 'number of Trials']
data = {c: [] for c in columns}

for cat, catExps in expIDsByCat.iteritems():
    for expName in catExps:
        if expName in processedResultsDF.index and processedResultsDF.loc[expName, 'Result']:
            print("Doing {} of category {}".format(expName, cat))
            rda = RawDataAnalyser(expName=expName, dirpath=NIXPath)
            freqSecNames = rda.getContStats()
            for freq, secNames in freqSecNames.iteritems():
                    data['Category'].append(cat)
                    data['Experiment ID'].append(expName)
                    data['Frequency (Hz)'].append(freq)
                    data['number of Trials'].append(len(secNames))

rawDF = pd.DataFrame(data)
catFreqGrouped = rawDF.groupby(['Category', 'Frequency (Hz)'])
statsDF = pd.DataFrame()
statsDF['(number of Exps., number of Trials)'] = catFreqGrouped['number of Trials'].agg(lambda x: (len(x), sum(x)))
statsDFUnstacked = statsDF.unstack()
statsDFTexTable = statsDFUnstacked.to_latex(column_format= 'c' + 'c' * len(statsDFUnstacked.columns))

statsDF2Write = rawDF.set_index(["Category", "Experiment ID", "Frequency (Hz)"])
statsDF2Write.to_excel(contStimDetailedStatsFile)


outDoc = pylatex.Document(contStimsResultsFile, documentclass='standalone', page_numbers=False, indent=False)
outDoc.packages.append(pylatex.Package('booktabs'))
outDoc.append(pylatex.utils.NoEscape(statsDFTexTable))
outDoc.generate_tex()

