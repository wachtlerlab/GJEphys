import numpy as np
import quantities as qu
from GJEphys.NEOFuncs import getSpikeRateIn, simpleFloat, getSpikeAmps
from GJEphys.doubleExpFitting import doubleExpFun


def spontAct1Sec(resp, spikes, xBest, xData, expUnits=qu.Hz):

    spikeRate = getSpikeRateIn(spikes, intervalStart=-1 * qu.s, intervalEnd=0 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def spontAct3Sec(resp, spikes, xBest, xData, expUnits=qu.Hz):

    spikeRate = getSpikeRateIn(spikes, intervalStart=-3 * qu.s, intervalEnd=0 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def initSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):

    spikeRate = getSpikeRateIn(spikes, intervalStart=0 * qu.s, intervalEnd=75 * qu.ms)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate


def laterSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    spikeRate = getSpikeRateIn(spikes, intervalStart=75 * qu.ms, intervalEnd=1 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def totalSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    spikeRate = getSpikeRateIn(spikes, intervalStart=0 * qu.ms, intervalEnd=1 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def reboundSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    spikeRate = getSpikeRateIn(spikes, intervalStart=1 * qu.s, intervalEnd=1.05 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate

def afterReboundSpikeRate(resp, spikes, xBest, xData, expUnits=qu.Hz):
    spikeRate = getSpikeRateIn(spikes, intervalStart=1.05 * qu.s, intervalEnd=2.00 * qu.s)

    spikeRate *= simpleFloat(expUnits / qu.Hz)

    return spikeRate


def maxExciNormed(resp, spikes, xBest, xData, expUnits=qu.dimensionless):

    if len(spikes) == 0 or any(np.isnan(xBest)):
        return np.nan
    else:
        avgSpikeHeight = np.mean(getSpikeAmps(resp, spikes.times))

        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        exciSig = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)

        return exciSig.max() / avgSpikeHeight


def maxInhiNormed(resp, spikes, xBest, xData, expUnits=qu.dimensionless):

    if len(spikes) == 0 or any(np.isnan(xBest)):
        return np.nan
    else:
        avgSpikeHeight = np.mean(getSpikeAmps(resp, spikes.times)).magnitude

        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        inhiSig = doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)

        return inhiSig.max() / avgSpikeHeight

def offsetNormed(resp, spikes, xBest, xData, expUnits=qu.dimensionless):
    if len(spikes) == 0 or any(np.isnan(xBest)):
        return np.nan
    else:
        avgSpikeHeight = np.mean(getSpikeAmps(resp, spikes.times))
        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        return offset / avgSpikeHeight

def inhiReleaseAt(resp, spikes, xBest, xData, expUnits=qu.ms):

    if any(np.isnan(xBest)):
        return np.nan
    else:
        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        inhiRelease = onset2 + 3 * taud2
        inhiRelease *= simpleFloat(expUnits / qu.ms)
        return inhiRelease

def exciInhiRatio(resp, spikes, xBest, xData, expUnits=qu.dimensionless):
    if any(np.isnan(xBest)):
        return np.nan
    else:
        [A1, A2, onset1, onset2, taur1, taud1, taur2, taud2, offset] = xBest
        inhiSig = doubleExpFun(xData, A2, A2, onset2, 1 / taur2, 1 / taud2)
        exciSig = doubleExpFun(xData, A1, A1, onset1, 1 / taur1, 1 / taud1)

        return exciSig.max() / inhiSig.max()

def firstSpikeLatency(resp, spikes, xBest, xData, expUnits=qu.ms):

    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes):
        firstSpikeLat = initSpikes[0]
        firstSpikeLat.units = expUnits
        return float(firstSpikeLat)
    else:
        return np.nan

def secondSpikeBISI(resp, spikes, xBest, xData, expUnits=qu.ms):
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes) >= 2:
        toReturn = initSpikes[1] - initSpikes[0]
        toReturn.units = expUnits
        return float(toReturn)
    else:
        return np.nan


def thirdSpikeBISI(resp, spikes, xBest, xData, expUnits=qu.ms):
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes) >= 3:
        toReturn = initSpikes[2] - initSpikes[1]
        toReturn.units = expUnits
        return float(toReturn)
    else:
        return np.nan

def fourthSpikeBISI(resp, spikes, xBest, xData, expUnits=qu.ms):
    initSpikes = spikes[(spikes > 0) & (spikes < 75 * qu.ms)]
    initSpikes.sort()

    if len(initSpikes) >= 4:
        toReturn = initSpikes[3] - initSpikes[2]
        toReturn.units = expUnits
        return float(toReturn)
    else:
        return np.nan

