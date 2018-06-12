'''
Contains dictionaries mapping certain feature name abbreviations to their longer forms.
'''

import quantities as qu

"""
Mappings for meta data
"""
mdFN = {
    'expID':        'Experiment ID',
    'freq':         'Frequency in Hz',
    'laborState':   'Labor State',
    'trialName':    'TrialName',
    'trialStart':   'Start Time of Trial (s)'
}

"""
Mappings for response features
"""
fFN = {
    'OnsetE':       'Onset of Excitation in ms',
    'OnsetI':       'Onset of Inhibition in ms',
    'taurE':        'Rise time constant \n of Excitation in ms',
    'taurI':        'Rise time constant \n of Inhibition in ms',
    'taudE':        'Decay time constant \n of Excitation in ms',
    'taudI':        'Decay time constant \n of Inhibition in ms',
    'AE':           'Amplitude constant\n of Excitation (unitless)',
    'AI':           'Amplitude constant\n of Inhibition (unitless)',
    'offset':       'Offset (unitless)'
}

"""
Mappings for more response features
"""
newFFN = {
    'spontFR1': 'Spontaneous Activity\n Rate for 1s \nBefore trial (spikes per second)',
    'spontFR3': 'Spontaneous Activity\n Rate for 3s \nBefore trial (spikes per second)',
    'initFR': 'Spiking rate\n in [0, 75]ms\n (spikes per second)',
    'laterFR': 'Spiking rate\n in [75, 1000]ms\n (spikes per second)',
    'totalFR': 'Spiking rate\n in [0, 1000]ms\n (spikes per second)',
    'reboundFR': 'Spiking rate\n in [1000, 1050]ms\n (spikes per second)',
    'afterReboundFR': 'Spiking rate\n in [1050, 2000]ms\n (spikes per second)',
    'exciMaxNormed': 'Normalized peak \n Excitation (unitless)',
    'inhiMaxNormed': 'Normalized peak \n Inhibition (unitless)',
    'offsetNormed': 'Normalized offset (unitless)',
    'inhiReleaseAt': 'Time of release \n from Inhibition in ms',
    'exciInhiRatio': 'Peak Excitation\n to  Peak Inhibition (unitless)',
    'firstSpikeLatency': 'First Spike\nLatency (ms)',
    'secondSpikeBISI': 'First ISI (ms)',
    'thirdSpikeBISI': 'Second ISI (ms)',
    'fourthSpikeBISI': 'Third ISI (ms)',
}

"""
Mappings between feature abbreviations and their units.
"""
newFFNUnits = {
    'spontFR1': qu.Hz,
    'spontFR3': qu.Hz,
    'initFR': qu.Hz,
    'laterFR': qu.Hz,
    'totalFR': qu.Hz,
    'reboundFR': qu.Hz,
    'afterReboundFR': qu.Hz,
    'exciMaxNormed': qu.dimensionless,
    'inhiMaxNormed': qu.dimensionless,
    'offsetNormed': qu.dimensionless,
    'inhiReleaseAt': qu.ms,
    'exciInhiRatio': qu.dimensionless,
    'firstSpikeLatency': qu.ms,
    'secondSpikeBISI': qu.ms,
    'thirdSpikeBISI': qu.ms,
    'fourthSpikeBISI': qu.ms,
}

"""
Mappings between feature name abbreviations and name of functions in GJEphys.additionalEphysFuncs.py implementing them.
"""
extraParamFuncs = {
    'spontFR1': 'spontAct1Sec',
    'spontFR3': 'spontAct3Sec',
    'initFR': 'initSpikeRate',
    'laterFR': 'laterSpikeRate',
    'totalFR': 'totalSpikeRate',
    'reboundFR': 'reboundSpikeRate',
    'afterReboundFR': 'afterReboundSpikeRate',
    'exciMaxNormed': 'maxExciNormed',
    'inhiMaxNormed': 'maxInhiNormed',
    'offsetNormed': 'offsetNormed',
    'inhiReleaseAt': 'inhiReleaseAt',
    'exciInhiRatio': 'exciInhiRatio',
    'firstSpikeLatency': 'firstSpikeLatency',
    'secondSpikeBISI': 'secondSpikeBISI',
    'thirdSpikeBISI': 'thirdSpikeBISI',
    'fourthSpikeBISI': 'fourthSpikeBISI',
}

"""
Legacy
"""
paramIndex = {'AE': 0, 'AI': 1, 'OnsetE': 2, 'OnsetI': 3, 'taurE': 4, 'taudE': 5, 'taurI': 6, 'taudI': 7, 'offset': 8}

"""
A subset of newFFN
"""
spikeFRSpikeTimesFNs = {
    'spontFR3': 'Spontaneous Activity\n Rate for 3s \nBefore trial (spikes per second)',
    'initFR': 'Spiking rate\n in [0, 75]ms\n (spikes per second)',
    'laterFR': 'Spiking rate\n in [75, 1000]ms\n (spikes per second)',
    'totalFR': 'Spiking rate\n in [0, 1000]ms\n (spikes per second)',
    'reboundFR': 'Spiking rate\n in [1000, 1050]ms\n (spikes per second)',
    'firstSpikeLatency': 'First Spike\nLatency (ms)',
    'secondSpikeBISI': 'First ISI (ms)',
    'thirdSpikeBISI': 'Second ISI (ms)',
    'fourthSpikeBISI': 'Third ISI (ms)',
}

"""
A subset of newFFNUnits
"""
spikeFRSpikeTimesFuncs = {
    'spontFR3': 'spontAct3Sec',
    'initFR': 'initSpikeRate',
    'laterFR': 'laterSpikeRate',
    'totalFR': 'totalSpikeRate',
    'reboundFR': 'reboundSpikeRate',
    'firstSpikeLatency': 'firstSpikeLatency',
    'secondSpikeBISI': 'secondSpikeBISI',
    'thirdSpikeBISI': 'thirdSpikeBISI',
    'fourthSpikeBISI': 'fourthSpikeBISI',
}

def getXBestFromCurrentData(currentData):
    '''
    Legacy function from fitting double exponential to response baseline.
    :param currentData:
    :return:
    '''
    xBest = [0 for x in range(len(paramIndex))]

    for k, v in paramIndex.iteritems():
        xBest[v] = currentData[fFN[k]]

    return xBest