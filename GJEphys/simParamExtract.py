import numpy as np
import os
from GJEMS.ephys.GNodeUpoadHelpers import *
from neo import Spike2IO
import quantities as qu
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from neo import Block, Segment,Event, NeoHdf5IO, Epoch

#***********************************************************************************************************************

def fitLine(x, y):

    assert len(x) == len(y)

    matA = np.ones([len(x), 2])
    matB = np.asarray(y).T

    matA[:, 0] = np.asarray(x).T

    m, c = np.linalg.lstsq(matA, matB)[0]

    return m, c

#***********************************************************************************************************************


class SimulationParameterExtracter():

    def __init__(self):
        fitDetails = dict(pOpts=None, bestWindow=None, timeConstants=None, errors=None)
        #start and stop times of steady state response
        self.tSS_Start = 200 * qu.ms
        self.tSS_Stop = 250 *qu.ms


    #*******************************************************************************************************************

    def parseSpike2Data(self, ephysFile, startStopTimes=None):

        assert  os.path.isfile(ephysFile), ephysFile + ' does not exist'
        self.ephysFile = ephysFile
        self.expName = os.path.split(ephysFile)[1].strip('.smr')
        spike2Reader = Spike2IO(self.ephysFile)
        dataBlock = spike2Reader.read()[0]

        entireVoltageSignal = dataBlock.segments[0].analogsignals[0]
        entireCurrentSignal = dataBlock.segments[0].analogsignals[2]

        self.currentMarkerTimes = dataBlock.segments[0].eventarrays[1].times
        currentStr = dataBlock.segments[0].eventarrays[1].annotations['extra_labels']
        self.currentsMarked = np.asarray([float(s.strip('nA')) for s in currentStr]) * 0.1 * qu.nA

        if startStopTimes is None:

            expStartTime = self.currentMarkerTimes[0] \
                           - 0.5 * (self.currentMarkerTimes[1] - self.currentMarkerTimes[0])
            expEndTime = self.currentMarkerTimes[-1] \
                         + 1.5 * (self.currentMarkerTimes[-1] - self.currentMarkerTimes[-2])
        else:
            expStartTime = startStopTimes[0] * qu.s
            expEndTime = startStopTimes[1] * qu.s

        expStartTime = max(expStartTime, entireVoltageSignal.t_start)
        expEndTime = min(expEndTime, entireVoltageSignal.t_stop)

        expStartIndex = int((expStartTime - entireVoltageSignal.t_start) * entireVoltageSignal.sampling_rate)
        expEndIndex = int((expEndTime - entireVoltageSignal.t_start) * entireVoltageSignal.sampling_rate)

        self.voltageSignal = entireVoltageSignal[expStartIndex:expEndIndex + 1] * 20 * qu.mV
        self.currentSignal = entireCurrentSignal[expStartIndex:expEndIndex + 1] * qu.mA

    #*******************************************************************************************************************

    def extractResponses(self):

        thresholdedSignal = np.asarray(self.currentSignal > 0.05, int)
        thresholdedSigDiff = np.diff(thresholdedSignal)

        positiveEdges = np.where(thresholdedSigDiff == 1)[0]
        negativeEdges = np.where(thresholdedSigDiff == -1)[0]

        if negativeEdges[0] < positiveEdges[0]:
            negativeEdges = negativeEdges[1:]

        if positiveEdges[len(positiveEdges) - 1] > negativeEdges[len(negativeEdges) - 1]:
            positiveEdges = positiveEdges[:len(positiveEdges) - 1]

        assert len(positiveEdges) == len(negativeEdges), 'Lengths of positive and negative edges don\'t match'

        self.voltageTraces = []
        self.currentAmps = []
        self.restingMembranePotentials = []

        for (start, stop) in zip(positiveEdges, negativeEdges):

            noSamples50ms = (50 * qu.ms * self.voltageSignal.sampling_rate).magnitude
            startToUse = max(start - noSamples50ms, 0)
            self.restingMembranePotentials.append(np.mean(self.voltageSignal[startToUse:start]))
            self.voltageTraces.append(self.voltageSignal[start: stop + 1])
            startTime = self.currentSignal.t_start + start * self.currentSignal.sampling_period
            # import ipdb
            # ipdb.set_trace()

            self.currentAmps.append(self.currentsMarked[max(np.where(self.currentMarkerTimes < startTime)[0])])

    #*******************************************************************************************************************

    def vizIndividualResponses(self):

        plt.figure()
        plt.show(block=False)
        for seg in self.CCData.segments:
            for (vtrace, iAmp) in zip(self.voltageTraces, self.currentAmps):
                plt.cla()
                plt.plot(vtrace)
                plt.title(iAmp)
                plt.draw()
                ch = 'n'
                ch = raw_input('Next(n)/Quit(q):')
                if ch == 'q':
                    break
        plt.close()

    #*******************************************************************************************************************

    def vizSingleAmpResp(self):

        cols = matplotlib.rcParams['axes.color_cycle']
        plt.figure()
        plt.show(block=False)
        plots = []
        amps = []
        for segInd in range(len(self.CCData.segments)):
            # plt.cla()

            seg = self.CCData.segments[segInd]
            for sigInd in range(len(seg.analogsignals)):
                vTrace = seg.analogsignals[sigInd]
                tVec = range(len(vTrace)) * vTrace.sampling_period
                iAmp = unicode2quantities(seg.events[sigInd].label)
                p, = plt.plot(tVec, vTrace,color=cols[segInd%len(cols)])
            plots.append(p)
            amps.append(iAmp)
            plt.legend(plots, amps)
            plt.draw()

            ch = 'n'
            ch = raw_input('Next(n)/Quit(q):')
            if ch == 'q':
                break

    #*******************************************************************************************************************

    def fitExp(self, xData, yData):

        funcToFit = lambda x, offset, amp, iTau: offset + amp * np.exp(-x * iTau)
        pOpt, pCov = curve_fit(funcToFit, xdata=xData, ydata=yData)
        return pOpt, pCov

    #*******************************************************************************************************************

    # def getTimeConstants(self):

        # funcToFit = lambda x, offset, amp, iTau: offset + amp * np.exp(-x * iTau)
        # errFunc = lambda p, x, y: sum((y - funcToFit(x, p[0], p[1], p[2])) ** 2) / float(len(y))

        # self.fitDetails['timeConstants'] = []
        # self.fitDetails['pOpts'] = []
        # self.fitDetails['bestWindow'] = []
        # self.fitDetails['errors'] = []
        #
        # stopIndices = np.arange(100, 1000, 100)
        # startIndex = 10
        #
        # for vtrace in self.voltageTraces:
        #
        #     print 'Fitting Exponential for trace' + str(self.voltageTraces.index(vtrace))
        #     pOpts = []
        #     errors = []
        #     if len(vtrace) > (min(stopIndices) - 10):
        #
        #         for stopIndex in stopIndices:
        #             yData = vtrace[startIndex:stopIndex + 1]
        #             xData = np.arange(len(yData)) * yData.sampling_period.magnitude
        #             yData = yData.magnitude
        #             pStart = [np.mean(yData), max(abs(yData)) - abs(yData[0]),
        #                       5.0 / (stopIndex - startIndex) / vtrace.sampling_period.magnitude]
        #
        #             pOpt, pCov = curve_fit(funcToFit, xdata=xData, ydata=yData, p0=pStart, maxfev=int(1e6))
        #             pOpt = pOpt.tolist()
        #             pOpt.extend([startIndex, stopIndex])
        #             errors.append(errFunc(pOpt, xData, yData))
        #             pOpts.append(pOpt)
        #
        #         # lowestErrorIndex = np.argmin(errors)
        #         lowestErrorIndex = 0
        #
        #         self.fitDetails['bestWindow'].append(lowestErrorIndex)
        #         self.fitDetails['timeConstants'].append(1 / pOpts[lowestErrorIndex][2])
        #         self.fitDetails['errors'].append(errors)
        #
        #     else:
        #         self.fitDetails['timeConstants'].append(None)
        #         self.fitDetails['bestWindow'].append(None)
        #
        #     self.fitDetails['pOpts'].append(pOpts)

    #*******************************************************************************************************************

    def getTimeConstants(self):

        funcToFit = lambda x,amp1,Itau1 : amp1 * (1 - np.exp(-x * Itau1))

        errFunc = lambda p, x, y : sum((y - funcToFit(x, p[0], p[1])) ** 2) / float(len(y))

        self.fitParams = []

        # plt.figure()
        # plt.show(block=False)

        for seg in self.CCData.segments:
            segFitParams = []
            for vTrace, epoch, event in zip(seg.analogsignals, seg.epochs, seg.events):
                restingPotential = unicode2quantities(epoch.label)
                startIndex = 10 + np.argmin(vTrace[10:100] < restingPotential)
                iAmp = unicode2quantities(event.label)
                yData = vTrace[startIndex:150] - vTrace[startIndex]
                xData = np.arange(len(yData)) * yData.sampling_period.magnitude
                yData = yData.magnitude
                pOpt, pCov = curve_fit(funcToFit, xdata=xData, ydata=np.sign(iAmp.magnitude) * yData, maxfev=int(1e6))
                segFitParams.append(pOpt)

                # plt.cla()
                # vTraceValid = vTrace[startIndex:]
                # allx = np.arange(len(vTraceValid)) * vTrace.sampling_period.magnitude
            #     plt.plot(allx, vTraceValid, 'r')
            #     plt.plot(allx, vTrace[startIndex].magnitude + np.sign(iAmp.magnitude) * np.array(funcToFit(allx,\
            #                                 pOpt[0], pOpt[1])), 'g')
            #     plt.title(str(1000/pOpt[1]) + 'ms')
            #     plt.draw()
            #
            #     ch = raw_input('Next(n)/Quit(q):')
            #     if ch == 'q':
            #         break
            #
            # if ch == 'q':
            #     break

            self.fitParams.append(segFitParams)

    #*******************************************************************************************************************

    def vizRinAndTaum(self):

        self.getTimeConstants()

        fig1 = plt.figure()
        plt.show(block=False)

        vDiffMeans = []
        iAmps = []
        timeConstants = [[1000/pOpt[1] for pOpt in segData] for segData in self.fitParams]

        for segInd in range(len(self.CCData.segments)):

            vDiffs = []
            seg = self.CCData.segments[segInd]

            for sigInd in range(len(seg.analogsignals)):
                vTrace = seg.analogsignals[sigInd]
                tVec = range(len(vTrace)) * vTrace.sampling_period
                iAmp = unicode2quantities(seg.events[sigInd].label)
                restingMemPot = unicode2quantities(seg.epochs[sigInd].label)
                indexStartSS = np.floor((self.tSS_Start * vTrace.sampling_rate).simplified.magnitude)
                indexEndSS = np.ceil((self.tSS_Stop * vTrace.sampling_rate).simplified.magnitude)
                voltageSS = np.mean(vTrace[indexStartSS:indexEndSS])

                if sum(np.shape(vTrace[indexStartSS:indexEndSS])):

                    vDiffs.append(voltageSS - restingMemPot)

                else:
                    print('No voltage Signal in the time range ' + str([self.tSS_Start, self.tSS_Stop])
                            + 'for segment=' + str(segInd) + 'trial=' + str(sigInd))



            plt.subplot(211)
            plt.errorbar(iAmp.magnitude, np.mean(vDiffs), yerr=np.std(vDiffs), color='r', marker='o', mfc='r', ms=5)
            plt.draw()

            vDiffMeans.append(np.mean(vDiffs))
            iAmps.append(iAmp.magnitude)

        m, c = fitLine(x=iAmps, y=vDiffMeans)


        iAmps = np.asarray(iAmps)


        plt.subplot(211)
        plt.plot(iAmps, m * iAmps + c, 'b')
        plt.title('m = Rin = ' + str(round(m, 2)) + 'megaohms; c=' + str(c) + 'mV')
        plt.xlabel('Current injected in nA')
        plt.ylabel('Voltage change induced in mV')

        plt.axis('tight')
        plt.draw()

        plt.subplot(212)
        plt.errorbar(iAmps, [np.mean(x) for x in timeConstants], yerr=[np.std(x) for x in timeConstants], color='r', marker='o', mfc='r', ms=5)
        plt.xlabel('Current injected in nA')
        plt.ylabel('Time Constants in ms')
        plt.draw()

    #*******************************************************************************************************************

    def plotFittedExp(self):

        funcToFit = lambda x,amp1,Itau1 : amp1 * (1 - np.exp(-x * Itau1))

        errFunc = lambda p, x, y : sum((y - funcToFit(x, p[0], p[1])) ** 2) / float(len(y))

        plt.figure()
        plt.show(block=False)
        plt.xlabel('time(s)')
        plt.ylabel('Voltage(mV)')

        for seg in self.CCData.segments:
            segFitParams = []
            for vTrace, epoch, event in zip(seg.analogsignals, seg.epochs, seg.events):
                restingPotential = unicode2quantities(epoch.label)
                startIndex = 10 + np.argmin(vTrace[10:100] < restingPotential)
                iAmp = unicode2quantities(event.label)
                yData = vTrace[startIndex:150] - vTrace[startIndex]
                xData = np.arange(len(yData)) * yData.sampling_period.magnitude
                yData = yData.magnitude
                pOpt, pCov = curve_fit(funcToFit, xdata=xData, ydata=np.sign(iAmp.magnitude) * yData, maxfev=int(1e6))
                segFitParams.append(pOpt)

                # plt.cla()
                vTraceValid = vTrace[startIndex:]
                allx = np.arange(len(vTraceValid)) * vTrace.sampling_period.magnitude
                plt.plot(allx, vTraceValid, 'r')
                plt.plot(allx, vTrace[startIndex].magnitude + np.sign(iAmp.magnitude) * np.array(funcToFit(allx,\
                                            pOpt[0], pOpt[1])), 'g')
                plt.title(event.label + '; ' + str(1000/pOpt[1]) + 'ms')

                plt.draw()

                ch = raw_input('Next(n)/Quit(q):')
                if ch == 'q':
                    break

            if ch == 'q':
                break



    #*******************************************************************************************************************

    def saveCCData(self):

        currentAmpsSet = np.sort(list(set([float(x.magnitude) for x in self.currentAmps]))).tolist()

        self.CCData = Block('Current Clamp Data')
        self.CCData.segments = [Segment(name='Current Of ' + unicode(iAmp) + 'nA') for iAmp in currentAmpsSet]
        for iAmp, vTrace in zip(self.currentAmps, self.voltageTraces):
            presSegInd = currentAmpsSet.index(iAmp)
            self.CCData.segments[presSegInd].analogsignals.append(vTrace)
            self.CCData.segments[presSegInd].events.append(Event(time=vTrace.t_start,label=unicode(iAmp)))
            self.CCData.segments[presSegInd].epochs.append(Epoch(time=vTrace.t_start - 50 * qu.ms,
                                                                 duration=50 * qu.ms,
                                                        label=unicode(self.restingMembranePotentials[presSegInd])))

        writer = NeoHdf5IO(os.path.join(os.path.split(self.ephysFile)[0], self.expName + '_CC.hdf5'))
        writer.write_block(self.CCData)
        writer.close()

    #*******************************************************************************************************************

    def doesCCDataExist(self):
        return self.CCData is not None

    #*******************************************************************************************************************

    def loadCCData(self, fileName):

        loader = NeoHdf5IO(fileName)
        self.expName = os.path.split(fileName)[1].strip('.smr')
        self.CCData = loader.read_block()
        loader.close()

    #*******************************************************************************************************************

#***********************************************************************************************************************

