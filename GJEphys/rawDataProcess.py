from neo import AnalogSignal
import matplotlib.pyplot as plt
from scipy.signal import argrelmin, argrelmax
import os
from GJEphys.neoNIXIO import addMultiTag, addTag, dataArray2AnalogSignal, property2qu, simpleFloat, \
                                 addQuantity2section, createPosDA
from GJEphys.stimParEstimation import estimateGJStimPars
import nixio as nix
import numpy as np
import quantities as qu
import logging

def getSpikesIn(spikeTimes, intervalStart, intervalEnd):

    return spikeTimes[(spikeTimes >= intervalStart) & (spikeTimes < intervalEnd)]

#***********************************************************************************************************************


class RawDataProcessor(object):

    #*******************************************************************************************************************

    def downSampleVibSignal(self):

        self.downSamplingFactor = int(round(self.vibrationSignal.sampling_rate / (4 * self.maximumFreq)))
        newSamplingRate = self.vibrationSignal.sampling_rate / self.downSamplingFactor
        downSamplingIndices = range(0, self.vibrationSignal.shape[0], self.downSamplingFactor)

        vibSigDownMag = self.vibrationSignal.magnitude[downSamplingIndices]
        vibSigDownMag -= np.median(vibSigDownMag)
        self.vibrationSignalDownStdDev = np.std(vibSigDownMag) * self.vibrationSignal.units

        temp = AnalogSignal(signal=vibSigDownMag,
                                                units=self.vibrationSignal.units,
                                                sampling_rate=newSamplingRate,
                                                t_start=self.vibrationSignal.t_start)
        self.vibrationSignalDown = temp.reshape((temp.shape[0], ))


    #*******************************************************************************************************************

    def downSampleVoltageSignal(self, downSamplingFactor=None):

        if downSamplingFactor is None:

            downSamplingFactor = self.downSamplingFactor

            newSamplingRate = self.voltageSignal.sampling_rate / downSamplingFactor
            downSamplingIndices = range(0, self.voltageSignal.shape[0], downSamplingFactor)

            voltageSignalDown = AnalogSignal(signal=self.voltageSignal.magnitude[downSamplingIndices],
                                                    units=self.voltageSignal.units,
                                                    sampling_rate=newSamplingRate,
                                                    t_start=self.voltageSignal.t_start)


            voltageSignalDown = voltageSignalDown.reshape((voltageSignalDown.shape[0],))
            return voltageSignalDown


    #*******************************************************************************************************************
    def __init__(self, expName, dirpath, readOnly=False, askShouldReprocess=True):


        self.maximumFreq = 700 * qu.Hz

        self.expName = expName

        self.spikesDetected = False

        if readOnly:
            self.nixFile = nix.File.open(os.path.join(dirpath, expName + '.h5'), nix.FileMode.ReadOnly)
        else:
            self.nixFile = nix.File.open(os.path.join(dirpath, expName + '.h5'), nix.FileMode.ReadWrite)

        rawDataBlock = self.nixFile.blocks['RawDataTraces']

        if not readOnly:

            for sec in self.nixFile.sections:

                if sec.name == 'VibrationStimulii-Processed':

                    if askShouldReprocess:
                        ch = raw_input('This file has already been processed. Reprocess?(y/n):')
                    else:
                        ch = 'y'

                    if ch != 'y':
                        exit('Aborted')

                    else:

                        for tagName in [tag.name for tag in rawDataBlock.tags]:
                            del rawDataBlock.tags[tagName]

                        for tagName in [tag.name for tag in rawDataBlock.multi_tags]:
                            del rawDataBlock.multi_tags[tagName]

                        for dataArrayName in [da.name for da in rawDataBlock.data_arrays]:
                            if dataArrayName not in ['VibrationStimulus', 'CurrentInput', 'MembranePotential']:
                                del rawDataBlock.data_arrays[dataArrayName]

                        del self.nixFile.sections['VibrationStimulii-Processed']


        self.vibrationSignal = dataArray2AnalogSignal(rawDataBlock.data_arrays['VibrationStimulus'])
        self.voltageSignal = dataArray2AnalogSignal(rawDataBlock.data_arrays['MembranePotential'])

    #*******************************************************************************************************************

    def plotVibEpoch(self, epochTimes, signal=None, points=False):

        # indStart = int(epochTimes[0] * qu.s * self.entireVibrationSignal.sampling_rate + self.recordingStartIndex)
        # indEnd = int(epochTimes[1] * qu.s * self.entireVibrationSignal.sampling_rate + self.recordingStartIndex)

        # epochTVec = self.entireVibrationSignal.t_start + np.arange(indStart, indEnd) * self.entireVibrationSignal.sampling_period

        # plt.plot(epochTVec, self.entireVibrationSignal[indStart:indEnd], 'g' + extra)

        self.downSampleVibSignal()

        indStart = int(epochTimes[0] * qu.s * self.vibrationSignalDown.sampling_rate)
        indEnd = int(epochTimes[1] * qu.s * self.vibrationSignalDown.sampling_rate)

        epochTVec = self.vibrationSignalDown.t_start + np.arange(indStart, indEnd) * self.vibrationSignalDown.sampling_period

        stimEnds = (np.array(self.stimEndInds)) / self.downSamplingFactor
        stimStarts = (np.array(self.stimStartInds)) / self.downSamplingFactor

        stimEndsPresent = [x * self.vibrationSignalDown.sampling_period + self.vibrationSignalDown.t_start
                           for x in stimEnds if indStart <= x <= indEnd]
        stimStartsPresent = [x * self.vibrationSignalDown.sampling_period + self.vibrationSignalDown.t_start
                             for x in stimStarts if indStart <= x <= indEnd]

        voltageSignalDown = self.downSampleVoltageSignal()

        extra = ''
        if points:
            extra = '*-'


        fig = plt.figure()
        plt.show(block=False)

        viblne, = plt.plot(epochTVec.magnitude, self.vibrationSignalDown[indStart:indEnd].magnitude, 'g' + extra)
        plt.xlabel('time (' + str(epochTVec.units) + ')')
        voltlne, = plt.plot(epochTVec.magnitude, voltageSignalDown[indStart:indEnd].magnitude, 'b' + extra)

        if len(stimEndsPresent):
            plt.stem(stimStartsPresent, 20 * np.ones(np.shape(stimStartsPresent)), 'k')
            plt.stem(stimEndsPresent, 20 * np.ones(np.shape(stimEndsPresent)), 'm')

        lnes = [viblne, voltlne]
        labels = ['vibration', 'membrane voltage']

        if signal is not None:

            newSamplingRate = signal.sampling_rate / self.downSamplingFactor
            downSamplingIndices = range(0, signal.shape[0], self.downSamplingFactor)

            signalDown = AnalogSignal(signal=signal.magnitude[downSamplingIndices],
                                                    units=signal.units,
                                                    sampling_rate=newSamplingRate,
                                                    t_start=signal.t_start)

            slne, = plt.plot(epochTVec, signalDown[indStart:indEnd], 'r' + extra)

            lnes.append(slne)
            labels.append('external signal')

        plt.plot(epochTVec.magnitude, 0.2 * np.ones(epochTVec.shape), 'y')

        plt.legend(lnes, labels)

        plt.draw()

    #*******************************************************************************************************************

    def getStimulusEpochs(self, minDur=5 * qu.ms, minGap=14 * qu.ms, vibThresholdSig=0.75 * qu.um):
        """
        Using direct thesholding of signal amplitude and adaptive threshold to determine a valid gap between stimuli.
        :return:
        """

        vibThresholdSig.units = self.vibrationSignalDown.units
        sigM = np.abs(self.vibrationSignalDown.magnitude)
        thresholdedSig = np.asarray(sigM > vibThresholdSig.magnitude, int)
        thresholdedSigDiff = np.diff(thresholdedSig)
        positiveEdges = np.where(thresholdedSigDiff == 1)[0]
        negativeEdges = np.where(thresholdedSigDiff == -1)[0]

        #making sure len(positiveEdges) == len(negativeEdges)

        # ignore a stimulus which started before the recording
        if negativeEdges[0] < positiveEdges[0]:
            negativeEdges = negativeEdges[1:]

        #ignore a stimulus which is still on at the end of the recording
        if positiveEdges[-1] > negativeEdges[-1]:
            positiveEdges = positiveEdges[:-1]

        minGapSamples = int((minGap / self.vibrationSignalDown.sampling_period).simplified)
        allEdges = [(n, p) for p, n in zip(positiveEdges[1:], negativeEdges[:- 1]) if p - n >= minGapSamples]

        self.stimStartInds = [positiveEdges[0]]
        self.stimEndInds = []

        for (end, start) in allEdges:
            self.stimStartInds.append(start)
            self.stimEndInds.append(end)

        self.stimEndInds.append(negativeEdges[len(negativeEdges) - 1])


        self.stimStartInds = np.array(self.stimStartInds) * self.downSamplingFactor
        self.stimEndInds = np.array(self.stimEndInds) * self.downSamplingFactor

        self.stimDur = (self.stimEndInds - self.stimStartInds) * self.vibrationSignal.sampling_period.simplified

        mask = (self.stimDur >= minDur)

        self.stimDur = self.stimDur[mask]
        self.stimStartInds = self.stimStartInds[mask]
        self.stimEndInds = self.stimEndInds[mask]

        self.stimStartTimes = self.vibrationSignal.t_start + \
                              self.stimStartInds * self.vibrationSignal.sampling_period.simplified
        self.stimEndTimes = self.vibrationSignal.t_start + \
                            self.stimEndInds * self.vibrationSignal.sampling_period.simplified

        self.stimInterval = np.zeros(np.shape(self.stimDur))
        self.stimInterval[:-1] = np.diff(self.stimStartInds) \
                                                        * self.vibrationSignal.sampling_period.simplified.magnitude
        self.stimInterval = qu.Quantity(self.stimInterval, self.vibrationSignal.sampling_period.simplified.units)

    #*******************************************************************************************************************

    def extractResponse(self):

        if not self.spikesDetected:

            self.detectSpikes()

        stimStarts = (self.stimStartInds) / self.downSamplingFactor
        stimEnds = (self.stimEndInds) / self.downSamplingFactor
        samplingRateDown = self.vibrationSignalDown.sampling_rate.simplified.magnitude

        self.stimAmps = []
        self.stimFreqs = []
        self.responseVTraces = []
        self.stimTraces = []
        self.stimDurs = []

        for (stD, endD, st, end, stT, endT) in \
                zip(stimStarts, stimEnds, self.stimStartInds, self.stimEndInds, self.stimStartTimes, self.stimEndTimes):

            stimDown = self.vibrationSignalDown[stD:endD + 1]
            stimDownFFT = np.fft.rfft(stimDown, n=2048)
            self.stimFreqs.append(np.argmax(np.abs(stimDownFFT)) * samplingRateDown / 2 / len(stimDownFFT))

            stimAS = self.vibrationSignal[st:end + 1]
            stim = stimAS.magnitude
            allAmps = stim[np.concatenate((argrelmin(stim)[0], argrelmax(stim)[0]))]

            self.stimAmps.append(np.abs(allAmps).mean())

            self.responseVTraces.append(self.voltageSignal[st:end + 1])

            self.stimTraces.append(stimAS - np.mean(stimAS))

            self.stimDurs.append((end - st + 1) * self.vibrationSignal.sampling_period.simplified.magnitude)

        self.stimFreqs = qu.Quantity(self.stimFreqs, self.vibrationSignalDown.sampling_rate.simplified.units)
        self.stimAmps = qu.Quantity(self.stimAmps, self.vibrationSignal.units)
        self.stimDurs = qu.Quantity(self.stimDurs, self.vibrationSignal.sampling_period.simplified.units)



    #*******************************************************************************************************************

    def showcaseSegments(self, start=0):

        fig = plt.figure()
        plt.show(block=False)


        for ind in range(start, len(self.stimTraces)):

            plt.cla()

            stim = self.stimTraces[ind]
            resp = self.responseVTraces[ind]
            freq = self.stimFreqs[ind]
            amp = self.stimAmps[ind]
            spikeTimes = getSpikesIn(self.spikeTimes, resp.t_start, resp.t_stop)


            tVec = resp.times
            plt.plot(tVec, stim, 'b')
            plt.plot(tVec, resp, 'r')
            plt.plot(spikeTimes, [max(resp)] * len(spikeTimes), 'g^')
            inter = (np.array([self.stimStartInds[ind], self.stimEndInds[ind]])) * self.vibrationSignal.sampling_period.simplified
            plt.title('Freq=' + str(np.round(freq, 2)) + ',  Amp=' + str(np.round(amp, 2)) + ',  Interval=' + str(np.round(inter, 2)))
            # plt.axis('off')
            plt.draw()

            ch = raw_input('Next(n) or Quit(q):')

            if ch == 'q':

                return



    #*******************************************************************************************************************

    def sortSegments(self):

        rawSec = self.nixFile.sections['VibrationStimulii-Raw']
        writtenFreqs = property2qu(rawSec.sections['ContinuousStimulii'].props['FrequenciesUsed'])
        self.hasPulses = False

        if any([x.name == 'PulseStimulii' for x in rawSec.sections]):
            writtenPulseDurs = property2qu(rawSec.sections['PulseStimulii'].props['PulseDurations'])
            writtenPulseIntervals = property2qu(rawSec.sections['PulseStimulii'].props['PulseIntervals'])
            self.pulseInput = {}
            self.hasPulses = True

        def nearestQu(a, a0):

            return a[np.abs(a-a0).argmin()]

        def nearestLowerEqQu(a, a0):

            temp = a[a <= a0]
            if temp.shape[0]:
                return nearestQu(temp, a0)
            else:
                return nearestQu(a, a0)

        self.contInput = {}

        self.backwardISIs = np.empty_like(self.stimDurs)
        self.forwardISIs = np.empty_like(self.stimDurs)

        self.backwardISIs[0] = 1 * self.stimDurs.units
        self.forwardISIs[-1] = 1 * self.stimDurs.units

        temp = np.diff(self.stimStartInds) * self.vibrationSignal.sampling_period.simplified
        self.backwardISIs[1:] = temp
        self.forwardISIs[:-1] = temp

        for stimInd in range(len(self.stimDurs)):
            stimTime = simpleFloat(self.stimStartInds[stimInd] * self.voltageSignal.sampling_period)

            #if it's a pulse stimulus
            if (self.stimDurs[stimInd] < 0.5 * qu.s) and self.hasPulses:

                nearestPulseDur = nearestLowerEqQu(writtenPulseDurs, self.stimDurs[stimInd])
                nearestPulseInterval = nearestQu(writtenPulseIntervals, self.forwardISIs[stimInd])

                key = (simpleFloat(nearestPulseDur), simpleFloat(nearestPulseInterval))

                # if it's the first pulse of the pulse train
                if (self.backwardISIs[stimInd] > 0.5 * qu.s) and (self.forwardISIs[stimInd] < 0.5 * qu.s):

                    diff = (nearestPulseInterval - self.forwardISIs[stimInd])
                    if diff > 6 * qu.ms:
                        print(
                            'Time: {}; Large difference between estimated and '
                            'written values of pulse interval'.format(stimTime)
                        )
                        print('Estimated PI: {}; Closest Written PI: {}; Diff: {}'.format(
                            self.forwardISIs[stimInd],
                            nearestPulseInterval, diff))
                        print('Freq: {}; Duration: {}; Amplitude: {}'.format(self.stimFreqs[stimInd],
                                                                             self.stimDurs[stimInd],
                                                                             self.stimAmps[stimInd]))
                        print('This response will be ignored')
                    else:
                        self.pulseInput.setdefault(key, []).append([stimInd])


                if key in self.pulseInput:
                    # for the intermediate pulses of the pulse train
                    if (self.backwardISIs[stimInd] < 0.5 * qu.s) and (self.forwardISIs[stimInd] < 0.5 * qu.s):



                        diff = (nearestPulseInterval - self.forwardISIs[stimInd])
                        if diff > 6 * qu.ms:
                            print(
                                'Time: {}; Large difference between estimated and '
                                'written values of pulse interval'.format(stimTime)
                            )
                            print('Estimated PI: {}; Closest Written PI: {}; Diff: {}'.format(
                                self.forwardISIs[stimInd],
                                nearestPulseInterval, diff))
                            print('Freq: {}; Duration: {}; Amplitude: {}'.format(self.stimFreqs[stimInd],
                                                                                 self.stimDurs[stimInd],
                                                                                 self.stimAmps[stimInd]))
                            print('This response will be ignored')
                        else:
                            self.pulseInput[key][-1].append(stimInd)

                    # for the last pulse of the pulse train
                    if (self.backwardISIs[stimInd] < 0.5 * qu.s) and (self.forwardISIs[stimInd] > 0.5 * qu.s):

                        nearestPulseInterval = nearestQu(writtenPulseIntervals, self.backwardISIs[stimInd])
                        key = (simpleFloat(nearestPulseDur), simpleFloat(nearestPulseInterval))

                        diff = (nearestPulseInterval - self.backwardISIs[stimInd])
                        if diff > 6 * qu.ms:
                            print('Time: {}; Large difference between estimated and '
                            'written values of pulse interval'.format(stimTime)
                                )
                            print('Estimated PI: {}; Closest Written PI: {}; Diff: {}'.format(
                                self.backwardISIs[stimInd],
                                nearestPulseInterval, diff))
                            print('Freq: {}; Duration: {}; Amplitude: {}'.format(self.stimFreqs[stimInd],
                                                                                 self.stimDurs[stimInd],
                                                                                 self.stimAmps[stimInd]))
                            print('This response will be ignored')
                        else:
                            self.pulseInput[key][-1].append(stimInd)


                # isolated pulse?
                else:
                    if (self.backwardISIs[stimInd] > 0.5 * qu.s) and (self.forwardISIs[stimInd] > 0.5 * qu.s):

                        print('FYI: Isolated Pulse at '
                                        + str(self.stimStartInds[stimInd] * self.voltageSignal.sampling_period))
                        print('Freq: {}; Duration: {}; Amplitude: {}'.format(self.stimFreqs[stimInd],
                                                                             self.stimDurs[stimInd],
                                                                             self.stimAmps[stimInd]))
                        print('This response will be ignored')



            # if it's a continuous stimulus
            else:

                nearestFreq = simpleFloat(nearestQu(writtenFreqs, self.stimFreqs[stimInd]))
                freqDiff = (nearestFreq * qu.Hz) - self.stimFreqs[stimInd]
                if freqDiff > 20 * qu.Hz:
                    print('Time: {}; Large difference between estimated and written frequencies'.format(stimTime))
                    print('Estimated Frequency: {}; Closest Written Freq: {}; Diff: {}'.format(self.stimFreqs[stimInd],
                                                                                               nearestFreq, freqDiff))
                    print('Freq: {}; Duration: {}; Amplitude: {}'.format(self.stimFreqs[stimInd],
                                                                         self.stimDurs[stimInd],
                                                                         self.stimAmps[stimInd]))
                    print('This response will be ignored')
                else:
                    self.contInput.setdefault(nearestFreq, []).append(stimInd)

        if self.hasPulses:
            for k, v in self.pulseInput.items():

                onlyTrain = len(v) == 1
                for ind, stims in enumerate(v):

                    if len(stims) == 1:

                        stimInd = stims[0]
                        print('FYI: Isolated Pulse at '
                              + str(self.stimStartInds[stimInd] * self.voltageSignal.sampling_period))
                        print('Freq: {}; Duration: {}; Amplitude: {}'.format(self.stimFreqs[stimInd],
                                                                             self.stimDurs[stimInd],
                                                                             self.stimAmps[stimInd]))
                        print('This response will be ignored')

                        self.pulseInput[k].pop(ind)
                        if onlyTrain:
                            self.pulseInput.pop(k)

    #*******************************************************************************************************************

    def detectSpikes(self):

        from .spikeDetection import SpikeDetector

        sd = SpikeDetector(self.voltageSignal)
        sd.filterButterworth()

        sd.thresholdAndDetect(0.5, 5)

        self.spikeTimes = sd.spikeTimes

        self.spikesDetected = True

    #*******************************************************************************************************************

    def addStimAndSpikes(self, tagCount, typ, start, extent, blk, refs, sec):


        addTag(
                    name='Tag' + str(tagCount + 1),
                    type=typ,
                    position=start,
                    extent=extent,
                    blk=blk,
                    refs=refs,
                    metadata=sec)

        spikes = getSpikesIn(self.spikeTimes, start, start + extent)


        addMultiTag(
                        name='MultiTag' + str(tagCount + 1),
                        type=typ,
                        positions=spikes,
                        blk=blk,
                        refs=refs,
                        metadata=sec
                   )

        return tagCount + 1

    #*******************************************************************************************************************

    def addContTags(self, stimInd, blk, sec, refs, tagCount, Fs):

        duringStart = simpleFloat(self.stimStartTimes[stimInd])
        beforeStart = simpleFloat(max(self.voltageSignal.t_start, self.stimStartTimes[stimInd] - Fs))
        afterStart = simpleFloat(self.stimEndTimes[stimInd])

        duringExtent = simpleFloat(self.stimDurs[stimInd])
        beforeExtent = duringStart - beforeStart
        afterExtent = simpleFloat(min(Fs, self.vibrationSignal.t_stop - self.stimEndTimes[stimInd]))

        tagCount = self.addStimAndSpikes(tagCount, 'DuringStimulus', duringStart, duringExtent, blk, refs, sec)
        tagCount = self.addStimAndSpikes(tagCount, 'BeforeStimulus', beforeStart, beforeExtent, blk, refs, sec)
        tagCount = self.addStimAndSpikes(tagCount, 'AfterStimulus', afterStart, afterExtent, blk, refs, sec)

        return tagCount


    #*******************************************************************************************************************

    def addPulseTags(self, stimInds, blk, sec, refs, tagCount, Fs):

        beforeStart = simpleFloat(max(0 * qu.s, self.stimStartTimes[stimInds[0]] - Fs))
        duringStart = simpleFloat(self.stimStartTimes[stimInds[0]])
        onPulseStarts = [simpleFloat(self.stimStartTimes[stimInd]) for stimInd in stimInds]
        offPulseStarts = [simpleFloat(self.stimEndTimes[stimInd]) for stimInd in stimInds]
        afterStart = simpleFloat(self.stimEndTimes[stimInds[-1]]) + onPulseStarts[-1] - offPulseStarts[-2]

        beforeExtent = duringStart - beforeStart
        duringExtent = afterStart - duringStart
        onPulseExtents = [simpleFloat(self.stimDurs[stimInd]) for stimInd in stimInds]
        offPulseExtents = [o - f for o, f in zip(onPulseStarts[1:], offPulseStarts[:-1])]
        offPulseExtents.append(offPulseExtents[-1])
        afterExtent = simpleFloat(min(Fs, self.vibrationSignal.t_stop - self.stimEndTimes[stimInds[-1]]))


        tagCount = self.addStimAndSpikes(tagCount, 'BeforeStimulus', beforeStart, beforeExtent, blk, refs, sec)
        tagCount = self.addStimAndSpikes(tagCount, 'DuringStimulus', duringStart, duringExtent, blk, refs, sec)


        for pulseInd in range(len(stimInds)):

            tagCount = self.addStimAndSpikes(tagCount, 'OnPulse' + str(pulseInd),
                                             onPulseStarts[pulseInd], onPulseExtents[pulseInd], blk, refs, sec)

            tagCount = self.addStimAndSpikes(tagCount, 'OffPulse' + str(pulseInd),
                                             offPulseStarts[pulseInd], offPulseExtents[pulseInd], blk, refs, sec)


        tagCount = self.addStimAndSpikes(tagCount, 'AfterStimulus', afterStart, afterExtent, blk, refs, sec)

        return tagCount

    #*******************************************************************************************************************

    def write2nixFile(self):

        procSec = self.nixFile.create_section('VibrationStimulii-Processed', 'Analysis')

        contSec = procSec.create_section('ContinuousStimulii', 'Stimulii/Sine')

        dataBlock = self.nixFile.blocks['RawDataTraces']
        vibSig = dataBlock.data_arrays['VibrationStimulus']
        voltSig = dataBlock.data_arrays['MembranePotential']

        tagCount = 0


        for freq, stimInds in self.contInput.iteritems():

            sec = contSec.create_section('ContinuousStimulusAt' + str(freq), 'Stimulus/Sine')
            addQuantity2section(sec, freq * qu.Hz, 'Frequency')

            for trialInd, stimInd in enumerate(stimInds):

                trialSec = sec.create_section('Trial' + str(trialInd + 1), 'StimulusTrial/Sine')
                addQuantity2section(trialSec, self.stimStartTimes[stimInd], 'StimulusStart')
                addQuantity2section(trialSec, self.stimEndTimes[stimInd], 'StimulusEnd')
                addQuantity2section(trialSec, self.stimAmps[stimInd], 'StimulusAmplitude')
                # addQuantity2section(trialSec, self.stimPhaseOffsets[stimInd], 'StimulusPhaseOffset')

                tagCount = self.addContTags(stimInd, dataBlock,
                                                    trialSec, [vibSig, voltSig], tagCount, 3 * qu.s)

        if self.hasPulses:

            pulseSec = procSec.create_section('PulseStimulii', 'Stimulii/Pulse')
            pulseInd = 1

            for pulsePar in self.pulseInput.iterkeys():

                sec = pulseSec.create_section('PulseStimulus' + str(pulseInd), 'Stimulus/PulseTrain')
                pulseInd += 1

                addQuantity2section(sec, 265 * qu.Hz, 'Frequency')
                addQuantity2section(sec, pulsePar[0] * self.stimDurs.simplified.units, 'PulseDuration')
                addQuantity2section(sec, pulsePar[1] * self.stimDurs.simplified.units, 'PulseInterval')


                for trialInd, stimInds in enumerate(self.pulseInput[pulsePar]):

                    trialSec = sec.create_section('Trial' + str(trialInd + 1), 'StimulusTrial/PulseTrain')
                    addQuantity2section(trialSec, self.stimStartTimes[stimInds[0]], 'StimulusStart')
                    addQuantity2section(trialSec, self.stimEndTimes[stimInds[-1]], 'StimulusEnd')
                    addQuantity2section(trialSec, self.stimAmps[stimInds[0]], 'StimulusAmplitude')
                    addQuantity2section(trialSec, len(stimInds) * qu.dimensionless, 'NumberOfPulseRepetitions')

                    tagCount = self.addPulseTags(stimInds, dataBlock,
                                                 trialSec, [vibSig, voltSig], tagCount, 3 * qu.s)

    #*******************************************************************************************************************

    def close(self):
        if self.nixFile.is_open():
            self.nixFile.close()

#***********************************************************************************************************************


