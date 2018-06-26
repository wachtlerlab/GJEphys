'''
This file contains functions to extract scalar features from eletrophysiological responses of DL-Int-1. All of them
have the same set of following input arguments:
:param resp:  neo.analogsignal, containing a response trace starting 3s before stimulus application and ending 3s
after stimulus application. The times of this trace need to be adjusted so that time at stimulus application is 0.
:param spikes: neo.spiketrain, containing the spikes corresponding to the interval of "resp" above. The spike times need
to be adjusted so that stimulus application is at time 0.
:param xBest: not relevant, legacy from double exponential curve fitting of response baseline.
See GJEphys/fitDoubleExp.py.
:param xData: not relevant, legacy from double exponential curve fitting of response baseline.
See GJEphys/fitDoubleExp.py.
:param expUnits: quantities.unitquantity.UnitQuantity, expected units for the quantity (Eg: quantities.Hz)

Not all inputs are used by each function. Those that are not used are indicated in documentation strings.
'''

import numpy as np
import quantities as qu
from GJEphys.NEOFuncs import getSpikeRateIn, simpleFloat, getSpikeAmps
from GJEphys.doubleExpFitting import doubleExpFun


def spontAct1Sec(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during the 1s interval preceding stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantities.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=-1 * qu.s, intervalEnd=0 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def spontAct3Sec(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during the 3s interval preceding stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantities.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=-3 * qu.s, intervalEnd=0 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def initSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during interval [0, 75)ms of stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=0 * qu.s, intervalEnd=75 * qu.ms)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate


def laterSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during interval [75, 1000)ms of stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=75 * qu.ms, intervalEnd=1 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def totalSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during interval [0, 1000)ms of stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=0 * qu.ms, intervalEnd=1 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def reboundSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during interval [25, 100)ms following the end of stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=1.025 * qu.s, intervalEnd=1.1 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def afterReboundSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    '''
    Calculates the spike rate during interval [100, 1000)ms following the end of stimulus application
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    spikeRate = getSpikeRateIn(spikes, intervalStart=1.1 * qu.s, intervalEnd=2.00 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate


def maxExciNormed(resp, spikes, xBest, xData, expUnits=qu.dimensionless):
    '''
    Legacy function from double fitting of response baseline.
    :param resp: see file documentation above
    :param spikes: see file documentation above
    :param xBest: see file documentation above
    :param xData: see file documentation above
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    if len(spikes) == 0 or any(np.isnan(xBest)):
        return np.nan
    else:
        avgSpikeHeight = np.mean(getSpikeAmps(resp, spikes.times))

        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        exciSig = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)

        return exciSig.max() / avgSpikeHeight


def maxInhiNormed(resp, spikes, xBest, xData, expUnits=qu.dimensionless):
    '''
    Legacy function from double fitting of response baseline.
    :param resp: see file documentation above
    :param spikes: see file documentation above
    :param xBest: see file documentation above
    :param xData: see file documentation above
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    if len(spikes) == 0 or any(np.isnan(xBest)):
        return np.nan
    else:
        avgSpikeHeight = np.mean(getSpikeAmps(resp, spikes.times)).magnitude

        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        inhiSig = doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)

        return inhiSig.max() / avgSpikeHeight

def offsetNormed(resp, spikes, xBest, xData, expUnits=qu.dimensionless):
    '''
    Legacy function from double fitting of response baseline.
    :param resp: see file documentation above
    :param spikes: see file documentation above
    :param xBest: see file documentation above
    :param xData: see file documentation above
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    if len(spikes) == 0 or any(np.isnan(xBest)):
        return np.nan
    else:
        avgSpikeHeight = np.mean(getSpikeAmps(resp, spikes.times))
        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        return offset / avgSpikeHeight

def inhiReleaseAt(resp, spikes, xBest, xData, expUnits=qu.ms):
    '''
    Legacy function from double fitting of response baseline.
    :param resp: see file documentation above
    :param spikes: see file documentation above
    :param xBest: see file documentation above
    :param xData: see file documentation above
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    if any(np.isnan(xBest)):
        return np.nan
    else:
        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        inhiRelease = onset2 + 3 * taud2
        inhiRelease *= simpleFloat(expUnits / qu.ms)
        return inhiRelease

def exciInhiRatio(resp, spikes, xBest, xData, expUnits=qu.dimensionless):
    '''
    Legacy function from double fitting of response baseline.
    :param resp: see file documentation above
    :param spikes: see file documentation above
    :param xBest: see file documentation above
    :param xData: see file documentation above
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    if any(np.isnan(xBest)):
        return np.nan
    else:
        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        inhiSig = doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)
        exciSig = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)

        return exciSig.max() / inhiSig.max()

def firstSpikeLatency(resp, spikes, xBest, xData, expUnits=qu.ms):
    '''
    Calculates the time interval from stimulus application to the first spike of the response
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes):
        firstSpikeLat = initSpikes[0]
        firstSpikeLat.units = expUnits
        return float(firstSpikeLat)
    else:
        return np.nan

def secondSpikeBISI(resp, spikes, xBest, xData, expUnits=qu.ms):
    '''
    Calculates the time interval between the first and second spikes of the response
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes) >= 2:
        toReturn = initSpikes[1] - initSpikes[0]
        toReturn.units = expUnits
        return float(toReturn)
    else:
        return np.nan


def thirdSpikeBISI(resp, spikes, xBest, xData, expUnits=qu.ms):
    '''
    Calculates the time interval between the second and third spikes of the response
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes) >= 3:
        toReturn = initSpikes[2] - initSpikes[1]
        toReturn.units = expUnits
        return float(toReturn)
    else:
        return np.nan

def fourthSpikeBISI(resp, spikes, xBest, xData, expUnits=qu.ms):
    '''
    Calculates the time interval between the third and fourth spikes of the response
    :param resp: unused
    :param spikes: see file documentation above
    :param xBest: unused
    :param xData: unused
    :param expUnits: see file documentation above
    :return: quantites.Quantity
    '''
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes) >= 4:
        toReturn = initSpikes[3] - initSpikes[2]
        toReturn.units = expUnits
        return float(toReturn)
    else:
        return np.nan

