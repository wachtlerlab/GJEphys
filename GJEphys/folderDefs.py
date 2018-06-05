import os
from GJEphys.KKHAXLParsing import parseMetaDataFile, getExpIDsByCategory

homeFolder = "/home/aj/"

excel = os.path.join(homeFolder, 'DataAndResults/Ginjang-Metadata/neuron_database.xlsx')
excelSheet = 'Kai-san final report150803'
NIXPath = os.path.join(homeFolder, 'DataAndResults/GJEphys/NIXFiles/')
spike2Path = os.path.join(homeFolder, 'DataAndResults/GJEphys/Ginjang-Ephys/')
contStimsResultsFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'Results', 'contStimStats')
contStimDetailedStatsFile = os.path.join(homeFolder, "DataAndResults", "GJEphys", "Results", "contStimDetailedStats.xlsx")
pulseStimsSummaryStatsFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'Results', 'pulseStimStats')
pulseStimsDetailedStatsFile = os.path.join(homeFolder, 'DataAndResults', 'GJEphys', 'Results', 'pulseStimDetailedStats.xlsx')
importResultsFile = os.path.join(NIXPath, "ImportResults.xlsx")
processedResultsFile = os.path.join(NIXPath, "processingResults.xlsx")

mdDF = parseMetaDataFile(excel, excelSheet, spike2Path)

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


def getAllExpNames():
    expNames = map(str, mdDF.index)
    return expNames

def allExpIDsByCat():

    return getExpIDsByCategory(excel, excelSheet, categories, spike2Path)