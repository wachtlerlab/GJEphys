# Ajayrama Kumaraswamy, 2016
# Ginjang project, LMU
"""
Contains functions and classes useful for extracting stimulus and response traces as well as response spikes as
neo analogsignals and spiketrains.
"""

from GJEphys.neoNIXIO import tag2AnalogSignal, multiTag2SpikeTrain, property2qu, getTagPosExt
import nixio as nix
import os
import quantities as qu
from neo import SpikeTrain
import numpy as np
import pandas as pd
from neoNIXIO import simpleFloat

def getSecPulsePars(sec):
    '''
    Converts and returns Pulse duration, pulse interval and frequency of pulse stimuli contained in properties of
    the section "sec"
    :param sec: an NIX section
    :return: tuple of three quantities.Quantity objects, (<pulse duration>, <pulse interval>, <frequency>)
    '''
    pulseDur = property2qu(sec.props['PulseDuration'])
    pulseDur.units = qu.ms

    pulseInt = property2qu(sec.props['PulseInterval'])
    pulseInt.units = qu.ms

    pulseFreq = property2qu(sec.props['Frequency'])
    pulseFreq.units = qu.Hz

    return (float(pulseDur.magnitude), float(pulseInt.magnitude), float(pulseFreq.magnitude))


def _getSpikeTrainStartStop(tagType, tagMetadata, Fs):
    '''
    Internal function, don't use
    '''

    stimStart = property2qu(tagMetadata.props['StimulusStart'])[0]
    stimEnd = property2qu(tagMetadata.props['StimulusEnd'])[0]

    if tagType == 'DuringStimulus':

        tStart = stimStart
        tStop = stimEnd

    if tagType == 'BeforeStimulus':

        tStop = stimStart
        tStart = tStop - Fs

    if tagType == 'AfterStimulus':

        tStart = stimEnd
        tStop = tStart + Fs

    return tStart, tStop


def isPulseTag(tag):
    """
    Returns True if the type of <tag> is "OnPulse" or "OffPulse"
    :param tag: nix.Tag object
    :return: bool
    """
    return tag.type.startswith('OnPulse') or tag.type.startswith('OffPulse')


def splitPulseTagType(tag):
    """
    Parses the type of <tag> to determine and return the type of pulse tag and the number associated with it. It
    requires that the type of <tag> is formatted as "<tag type><pulse number>" where <tag type> is "OnPulse" or
    "OffPulse" and <pulse number> is an integer.
    :param tag: nix.Tag object
    :return: (int, string), (pulse number, pulse stimulus type)
    """

    assert isPulseTag(tag), 'Tag {} of type {} is not a pulse tag'.format(tag.name, tag.type)

    if tag.type.startswith('OnPulse'):

        pulseNumber = int(tag.type[7:])
        pulseStimType = 'OnPulses'

    elif tag.type.startswith('OffPulse'):
        pulseNumber = int(tag.type[8:])
        pulseStimType = 'OffPulses'

    else:
        raise(Exception('Tag {} of type {} is not a pulse tag'.format(tag.name, tag.type)))

    return pulseNumber, pulseStimType


def secPres(testSec, secs):
    """
    Legacy function, do not use.
    :param testSec:
    :param secs:
    :return:
    """

    for secInd, sec in enumerate(secs):
        freqPres = testSec.props['Frequency'].values[0].value == sec.props['Frequency'].values[0].value
        pulDurPres = testSec.props['PulseDuration'].values[0].value == sec.props['PulseDuration'].values[0].value
        pulIntPres = testSec.props['PulseInterval'].values[0].value == sec.props['PulseInterval'].values[0].value

        if freqPres and pulDurPres and pulIntPres:
            return secInd

    return -1



class RawDataAnalyser(object):
    '''
    Class used to collect methods for accessing the processed electrophysiology data of our project
    '''

    def __init__(self, expName, dirpath):
        '''
        :param expName: string, experiment name/ID
        :param dirpath: string, path for finding the NIX files
        '''

        self.expName = expName
        self.nixFile = nix.File.open(str(os.path.join(dirpath, expName + '.h5')), nix.FileMode.ReadOnly)

    def getContStats(self):
        """
        Returns a dictionary with frequencies used for continuous stimuli as keys and lists of associated trial names
        as values. Trial names are formatted as "TrialX" where X is the trial number of the associated stimulus
        frequency.
        :return: dict
        """

        contSec = self.nixFile.sections['VibrationStimulii-Processed'].sections['ContinuousStimulii']
        freqStats = {}
        for sec in contSec.sections:
            freq = sec.props['Frequency'].values[0].value
            freqStats[freq] = [s.name for s in sec.sections]

        return freqStats

    def getPulseStats(self):
        """
        Returns a dictionary with combinations of pulse parameters, formatted as (<Pulse Duration>, <Pulse Interval>),
        as keys and corresponding lists of trial names as values. Trial names are formatted as "TrialX" where X is the
        trial number of the associated pulse parameter combination.
        :return: dict
        """
        tempDF = pd.DataFrame()
        if 'PulseStimulii' in self.nixFile.sections['VibrationStimulii-Processed'].sections:
            pulseParSecNames = {}
            for sec in self.nixFile.sections['VibrationStimulii-Processed'].sections['PulseStimulii'].sections:
                pulsePars = getSecPulsePars(sec)
                for s in sec.sections:
                    tempS = pd.Series()
                    tempS['Pulse Duration (ms)'] = pulsePars[0]
                    tempS['Pulse Interval (ms)'] = pulsePars[1]
                    tempS['Trial Label'] = s.name
                    tempS['Frequency (Hz)'] = pulsePars[2]
                    tempS['Number of Pulse Repetitions'] = \
                        float(property2qu(s.props['NumberOfPulseRepetitions']).magnitude)
                    tStart = property2qu(s.props["StimulusStart"])[0]
                    tStart.unit = qu.s
                    tempS['Time of Start of Pulse train (s)'] = simpleFloat(tStart)
                    tempDF = tempDF.append(tempS, ignore_index=True)
        return tempDF


    def getContResps(self, freqs=None, types=None):
        '''
        Collects and returns the responses from the current experiment to vibration stimulii of specified frequencies and types
        :param freqs: iterable of floats, the frequencies of vibration stimulii to which responses are collected.
        Default: List of all frequencies used in the experiment
        :param types: list, must be a subset of ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus']. The strings
        correspond to time intervals of 3 seconds stimulus onset, during the stimulus and 3 seconds after stimulus
        offset. Default: ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus'].
        :return: dict, with three levels of nesting. First level of nesting has combination of pulse parameters
        as keys. The second level has trial names as keys. The third level has strings in <types> as keys and response
        traces formatted as neo.analogsignal as values.
        '''

        if freqs is None:

            sec = self.nixFile.sections['VibrationStimulii-Raw'].sections['ContinuousStimulii']
            freqs = [v.value for v in sec.props['FrequenciesUsed'].values]

        if types is None:

            types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']

        resps = {}
        allFreqSecs = self.nixFile.sections['VibrationStimulii-Processed'].sections['ContinuousStimulii'].sections
        freqSecs = {}

        for fs in allFreqSecs:

            freq = fs.props['Frequency'].values[0].value
            if freq in freqs:
                freqSecs[freq] = fs


        for freq in freqs:

            if freq in freqSecs.viewkeys():
                resps[freq] = []
                nTrials = len(freqSecs[freq].sections)

                for ind in range(nTrials):

                    resps[freq].append({})

        freqSecs = freqSecs.values()


        for tag in self.nixFile.blocks['RawDataTraces'].tags:

            if tag.metadata.parent in freqSecs and tag.type in types and tag.metadata.type == 'StimulusTrial/Sine':

                freq = tag.metadata.parent.props['Frequency'].values[0].value

                analogSignal = tag2AnalogSignal(tag, 1)

                resps[freq][int(tag.metadata.name[5:]) - 1][tag.type] = analogSignal

        return resps

    def getContSpikes(self, freqs=None, types=None):
        '''
        Collects and returns the spike times during the responses of the current experiment to vibration stimulii of specified frequencies and types
        :param freqs: iterable of floats, the frequencies of vibration stimulii to which responses are collected.
        Default: List of all frequencies used in the experiment
        :param types: list, must be a subset of ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus']. The strings
        correspond to time intervals of 3 seconds stimulus onset, during the stimulus and 3 seconds after stimulus
        offset. Default: ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus'].
        :return: dict, with three levels of nesting. First level of nesting has frequency values as keys.
        The second level has trial names as keys. The third level has strings in <types> as keys and spike times
         formatted as neo.spiketrain as values.
        '''

        if freqs is None:

            sec = self.nixFile.sections['VibrationStimulii-Raw'].sections['ContinuousStimulii']
            freqs = [v.value for v in sec.props['FrequenciesUsed'].values]

        if types is None:

            types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus']

        spikes = {}

        allFreqSecs = self.nixFile.sections['VibrationStimulii-Processed'].sections['ContinuousStimulii'].sections
        freqSecs = {}

        Fs = 3 * qu.s

        for fs in allFreqSecs:

            freq = fs.props['Frequency'].values[0].value
            if freq in freqs:
                freqSecs[freq] = fs
                spikes[freq] = []
                trialSecs = [s for s in fs.sections if s.type == 'StimulusTrial/Sine']
                for sec in trialSecs:
                    temp = {}
                    for typ in types:

                        tStart, tStop = _getSpikeTrainStartStop(typ, sec, Fs)
                        temp[typ] = SpikeTrain(times=[], t_start=tStart, t_stop=tStop, units=tStart.units)
                    spikes[freq].append(temp)


        freqSecs = freqSecs.values()

        for tag in self.nixFile.blocks['RawDataTraces'].multi_tags:

            if tag.metadata.parent in freqSecs and tag.type in types and tag.metadata.type == 'StimulusTrial/Sine':

                freq = tag.metadata.parent.props['Frequency'].values[0].value

                tStart, tStop = _getSpikeTrainStartStop(tag.type, tag.metadata, Fs)

                sp = multiTag2SpikeTrain(tag, tStart, tStop)

                spikes[freq][int(tag.metadata.name[5:]) - 1][tag.type] = sp



        return spikes



    def getPulseResps(self, types=None, pulsePars=None):
        '''
        Collects and returns the responses from the current experiment to pulse stimulii of specified
        pulse parameter combinations.
        :param pulsePars: iterable of 2 float tuples, the combination of pulse parameters, for which responses are
        collected and returned. The tuples need to be formatted as (<Pulse Duration>, <Pulse Interval>).
        Default: All pulse parameter combinations used for the experiment.
        :param types: list, must be a subset of ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus']. The strings
        correspond to time intervals of 3 seconds stimulus onset, during the stimulus and 3 seconds after stimulus
        offset. Default: ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus'].
        :return: dict, with three levels of nesting. First level of nesting has combination of pulse parameters
        as keys. The second level has trial names as keys. The third level has strings in <types> as keys and response
        traces formatted as neo.analogsignal as values.
        '''

        if not 'PulseStimulii' in self.nixFile.sections['VibrationStimulii-Processed'].sections:

            return {}
        else:
            ppSec = self.nixFile.sections['VibrationStimulii-Processed'].sections['PulseStimulii']

            if types is None:
                types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus', 'OnPulses', 'OffPulses']

            if pulsePars is None:

                pulsePars = [getSecPulsePars(s) for s in ppSec.sections]

            resps = {}

            for sec in ppSec.sections:

                pp = getSecPulsePars(sec)
                if pp in pulsePars:
                    resps[pp] = [{k: {} for k in types} for s in sec.sections if s.type == 'StimulusTrial/PulseTrain']


            for tag in self.nixFile.blocks['RawDataTraces'].tags:

                tagTypeValid = (tag.type in types) or (tag.type.find('Pulse') >= 0)

                if tag.metadata.parent in ppSec.sections \
                                        and tagTypeValid and tag.metadata.type == 'StimulusTrial/PulseTrain':

                    pp = getSecPulsePars(tag.metadata.parent)
                    trialInd = int(tag.metadata.name[5:]) - 1

                    if isPulseTag(tag):
                        pulseNumber, pulseStimType = splitPulseTagType(tag)
                        resps[pp][trialInd][pulseStimType][pulseNumber] = tag2AnalogSignal(tag, 1)
                    else:
                        resps[pp][trialInd][tag.type] = tag2AnalogSignal(tag, 1)



            for pp, ppTrials in resps.iteritems():

                for trialInd, ppTrial in enumerate(ppTrials):

                    for typ, typResp in ppTrial.iteritems():

                        if typ in ['OnPulses', 'OffPulses']:

                            temp = [0 for x in xrange(max(typResp.keys()) + 1)]
                            for k, v in typResp.iteritems():
                                temp[k] = v

                            resps[pp][trialInd][typ] = temp

            return resps

    def getPulseSpikes(self, types=None, pulsePars=None):
        '''
       Collects and returns the spike times during the responses of the current experiment to pulse stimulii of specified
       pulse parameter combinations.
       :param pulsePars: iterable of 2 float tuples, the combination of pulse parameters, for which responses are
       collected and returned. The tuples need to be formatted as (<Pulse Duration>, <Pulse Interval>).
       Default: All pulse parameter combinations used for the experiment.
       :param types: list, must be a subset of ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus']. The strings
       correspond to time intervals of 3 seconds stimulus onset, during the stimulus and 3 seconds after stimulus
       offset. Default: ['BeforeStimulus', 'DuringStimulus', 'AfterSimulus'].
       :return: dict, with three levels of nesting. First level of nesting has combination of pulse parameters
       as keys. The second level has trial names as keys. The third level has strings in <types> as keys and spike trains
       formatted as neo.spiketrains as values.
       '''

        if not 'PulseStimulii' in self.nixFile.sections['VibrationStimulii-Processed'].sections:

            return {}
        else:
            ppSec = self.nixFile.sections['VibrationStimulii-Processed'].sections['PulseStimulii']

            if types is None:
                types = ['BeforeStimulus', 'DuringStimulus', 'AfterStimulus', 'OnPulses', 'OffPulses']

            if pulsePars is None:

                pulsePars = [getSecPulsePars(s) for s in ppSec.sections]

            spikes = {}

            for sec in ppSec.sections:

                pp = getSecPulsePars(sec)
                if pp in pulsePars:
                    spikes[pp] = [{k: {} for k in types} for s in sec.sections if s.type == 'StimulusTrial/PulseTrain']


            for tag in self.nixFile.blocks['RawDataTraces'].tags:

                tagTypeValid = (tag.type in types) or (tag.type.find('Pulse') >= 0)

                if tag.metadata.parent in ppSec.sections \
                                        and tagTypeValid and tag.metadata.type == 'StimulusTrial/PulseTrain':

                    pp = getSecPulsePars(tag.metadata.parent)
                    trialInd = int(tag.metadata.name[5:]) - 1

                    pos, ext = getTagPosExt(tag)

                    expectedName = 'MultiTag' + tag.name[3:]
                    if expectedName in self.nixFile.blocks['RawDataTraces'].multi_tags:
                        corrTag = self.nixFile.blocks['RawDataTraces'].multi_tags[expectedName]
                        spikeTrain = multiTag2SpikeTrain(corrTag, tStart=pos, tStop=pos + ext)
                    else:
                        spikeTrain = SpikeTrain([], t_start=pos, t_stop=pos+ext, units=pos.units)

                    if isPulseTag(tag):
                        pulseNumber, pulseStimType = splitPulseTagType(tag)
                        spikes[pp][trialInd][pulseStimType][pulseNumber] = spikeTrain
                    else:
                        spikes[pp][trialInd][tag.type] = spikeTrain



            for pp, ppTrials in spikes.iteritems():

                for trialInd, ppTrial in enumerate(ppTrials):

                    for typ, typResp in ppTrial.iteritems():

                        if typ in ['OnPulses', 'OffPulses']:

                            temp = [0 for x in xrange(max(typResp.keys()) + 1)]
                            for k, v in typResp.iteritems():
                                temp[k] = v

                            spikes[pp][trialInd][typ] = temp

            return spikes

