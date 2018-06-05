'''
Description: The script generates a table summarizing the number of trials and number of experiments for every neuron
category and for every combination of (Pulse Interval, Pulse Duration) of the stimulus applied.
This script uses NIX files in the directory "NIXPath" in GJEphy.folderDefs. The original smr files should be imported
using "rawDataImporter.py" and the resulting files processed using "rawDataProcessor.py" before running this script
to generate a table of summary statistics. This script checks whether the NIX file of an input <Experiment ID> has been
successfully processed by looking at the excel file at the path indicated by GJEphys.folderDefs.processedResultsFile.

Usage: python <path to parent directory>/pulseStimStats.py
'''

from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory
from GJEphys.rawDataAnalyse import RawDataAnalyser
import pylatex
import pandas as pd
from GJEphys.folderDefs import excel, excelSheet, NIXPath, spike2Path, pulseStimsSummaryStatsFile, \
    processedResultsFile, pulseStimsDetailedStatsFile

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
columns = ['Category', 'Exp ID', '(duration, interval) (ms)', 'number of Trials']
rawDF = pd.DataFrame()

for cat, catExps in expIDsByCat.iteritems():
    for expName in catExps:
        if expName in processedResultsDF.index and processedResultsDF.loc[expName, 'Result']:
            print("Doing {} of category {}".format(expName, cat))
            rda = RawDataAnalyser(expName=expName, dirpath=NIXPath)
            expPulseStatsDF = rda.getPulseStats()
            expPulseStatsDF['Experiment ID'] = expName
            expPulseStatsDF['Category'] = cat
            rawDF = rawDF.append(expPulseStatsDF)


rawDF['(Pulse Interval, Pulse Duration) (ms)'] = [(x, y) for x, y in
                                                             zip(rawDF['Pulse Interval (ms)'],
                                                                 rawDF['Pulse Duration (ms)'])]
del rawDF['Pulse Duration (ms)']
del rawDF['Pulse Interval (ms)']
catFreqGrouped = rawDF.groupby(['(Pulse Interval, Pulse Duration) (ms)', 'Category'])
statsDF = pd.DataFrame()
statsDF['(number of Exps., number of Trials)'] = \
    catFreqGrouped['Experiment ID'].agg(lambda x: (len(x.unique()), len(x)))
statsDFUnstacked = statsDF.unstack()
statsDFTexTable = statsDFUnstacked.to_latex()

detailedStatsDF = rawDF.set_index(keys=["Category", "Experiment ID", "Trial Label"])
detailedStatsDF.to_excel(pulseStimsDetailedStatsFile)

outDoc = pylatex.Document(pulseStimsSummaryStatsFile, documentclass='standalone', page_numbers=False, indent=False)
outDoc.packages.append(pylatex.Package('booktabs'))
outDoc.append(pylatex.utils.NoEscape(statsDFTexTable))
outDoc.generate_tex()

