import os
import numpy as np
import pandas as pd
import quantities as qu
import logging

#***********************************************************************************************************************

def smrFilesPresent(expNames, smrDir):
    temp1 = []
    dirList = os.listdir(smrDir)
    for expName in expNames:
        matches = [x.endswith('.smr') and x.startswith(expName) for x in dirList]
        smrFilePresent = any(matches)
        temp1.append(smrFilePresent)
    temp1 = np.array(temp1).reshape((len(temp1), 1))
    return temp1

#***********************************************************************************************************************

def addDyeNamesIfNess(expNames, smrDir):
    newIndices = {}
    dirList = os.listdir(smrDir)
    for expName in expNames:
        matches = [x.endswith('.smr') and x.startswith(expName) for x in dirList]
        if any(matches):
            
            newIndex = dirList[matches.index(True)][:-4]
            if len(newIndex) - len(expName) == 2:
                newIndices[expName] = newIndex
            else:
                newIndices[expName] = expName

        if expName not in newIndices:
            newIndices[expName] = expName

    return newIndices

#***********************************************************************************************************************
def parseMetaDataFile(excelFile, excelSheet, smrDir):

    tempDf = pd.read_excel(excelFile, sheetname=excelSheet, header=None, parse_cols=None)

    currentVal = tempDf.loc[0, 0]
    for ind, val in enumerate(tempDf.loc[0, 1:]):
        if pd.isnull(val):
            tempDf.loc[0, ind + 1] = currentVal
        else:
            currentVal = val

    for ind, val in enumerate(tempDf.loc[1, :]):

        if pd.isnull(val):
            tempDf.loc[1, ind] = tempDf.loc[0, ind]

    metaDF = tempDf.loc[3:, 1:]
    metaDF.index = tempDf.loc[3:, 0]
    metaDF.columns = pd.MultiIndex.from_arrays((tempDf.loc[0, 1:], tempDf.loc[1, 1:]))
    metaDF.sort_index()

    smrsPresent = smrFilesPresent(metaDF.index, smrDir)

    problem = np.logical_not(smrsPresent)

    if any(problem):
        logging.warning('These experiments were marked uploaded but were not found in {}:\n{}\nIGNORING THESE'
                        .format(smrDir, [x for x, y in zip(metaDF.index, problem) if y]))

    metaDFFiltered = metaDF[smrsPresent]


    expNames = [str(x) for x in metaDFFiltered.index]
    newIndices = addDyeNamesIfNess(expNames, smrDir)
    metaDFFiltered = metaDFFiltered.rename(index=newIndices)


    return metaDFFiltered

#***********************************************************************************************************************

def extractMetaData(mdDF, expName):

    expData = mdDF.loc[expName]

    freqEntry = str(expData[('Stimulus', 'Frequency')])
    pulseEntry = expData[('Stimulus', 'Pulse (Duration/Interval)')]
    spontEntry = expData[('Activity', 'Spontaneous')]
    respEntry = expData[('Activity', 'Response')]
    recPeriod = expData[('Recording period (s)', 'Recording period (s)')]
    voltCalibStr = expData[('calibration', 'mV/V (uploaded files)')]
    stimCalibStr = expData[('calibration', 'um/V')]
    currCalibStr = expData[('calibration', 'nA/V')]
    int2Exclude = expData[("Intervals To Exclude (s)", "Intervals To Exclude (s)")]



    metaData = {}
    metaData['pulse'] = [[], []]
    metaData['freqs'] = []
    metaData['spont'] = ''
    metaData['resp'] = ''
    metaData['recPeriod'] = None
    metaData['voltCalibStr'] = None
    metaData['stimCalibStr'] = None
    metaData['currCalibStr'] = None
    metaData['int2Exclude'] = None

    if not (pd.isnull(recPeriod) or recPeriod == '-'):
        metaData['recPeriod'] = recPeriod

    if not (pd.isnull(voltCalibStr) or voltCalibStr == '-'):
        metaData['voltCalibStr'] = str(voltCalibStr)

    if not (pd.isnull(stimCalibStr) or stimCalibStr == '-'):
        metaData['stimCalibStr'] = str(stimCalibStr)

    if not (pd.isnull(currCalibStr) or currCalibStr == '-'):
        metaData['currCalibStr'] = str(currCalibStr)

    if not (pd.isnull(int2Exclude) or int2Exclude == '-'):
        metaData['int2Exclude'] = int2Exclude


    if freqEntry != 'nan':
        if freqEntry.find(',') < 0:

            # sometimes the commas are not read in, leading to concatenated entries like 100265. Splitting them.
            tempN = len(freqEntry)
            for x in xrange(int(np.ceil(tempN / 3.0))):

                metaData['freqs'].append(int(freqEntry[max(0, tempN - (x + 1) * 3): tempN - 3 * x]))

        else:
            metaData['freqs'] = map(lambda x: int(x), freqEntry.split(','))

    metaData['freqs'] *= qu.Hz

    if not pd.isnull(pulseEntry):
        unresolved = []
        # pulseEntry is expected to be made of two types of entries:
        # (i) a/b (ii) a, c / b which is the same as a/b , c/b
        try:
            for word in pulseEntry.split(','):
                if word.count('/'):
                    (duration, interval) = word.split('/')
                    unresolved.append(float(duration))
                    metaData['pulse'][1].extend([float(interval)] * len(unresolved))
                    metaData['pulse'][0].extend(unresolved)
                    unresolved = []
                else:
                    unresolved.append(float(word))

            metaData['pulse'][1].extend(unresolved)
            lastNum = metaData['pulse'][0][-1]
            metaData['pulse'][0].extend([lastNum] * len(unresolved))
        except:
            raise(Exception('Improper entry in pulse column for the given smr file.'))

    metaData['pulse'][0] *= qu.ms
    metaData['pulse'][1] *= qu.ms

    if not pd.isnull(spontEntry):
        metaData['spont'] = bool(spontEntry.count('yes') + spontEntry.count('Yes') + spontEntry.count('YES'))


    if not pd.isnull(respEntry):
        metaData['resp'] = str(respEntry)

    return metaData

# **********************************************************************************************************************

def validExpName(expName):

    lenCheckPass = len(expName) in [8, 10]
    if not lenCheckPass:
        return False
    digitCheckPass = expName[:6].isdigit() and expName[7].isdigit() and expName[6] == '-'
    if not digitCheckPass:
        return False

    if len(expName) == 10:
        dyeCheckPass = expName[8:] in ['LY', 'Rh', 'Al']

        if not dyeCheckPass:
            return False

    return True
# **********************************************************************************************************************

def getExpIDsByCategory(excel, originalSheet, categorySheets, smrDir):

    mdDF = parseMetaDataFile(excel, originalSheet, smrDir)
    expIDsByCat = {}

    for category in categorySheets:

        categoryExps = pd.read_excel(excel, category, parse_cols=[0], header=None)
        categoryExps = filter(validExpName, map(str, categoryExps.iloc[:, 0]))

        categoryExps = addDyeNamesIfNess(categoryExps, smrDir).values()
        categoryExpsWithData = filter(lambda x: x in mdDF.index, categoryExps)

        expIDsByCat[category] = categoryExpsWithData

    return expIDsByCat
