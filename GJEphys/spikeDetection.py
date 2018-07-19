from scipy.signal import iirfilter, freqz, lfilter, kaiserord, firwin
from neo import AnalogSignal
from matplotlib import pyplot as plt
from math import pi as PI
import numpy as np
import quantities as qu
from sklearn.mixture import GaussianMixture
from scipy.stats import norm


def getLocalMaximaThresholded(a, thresh):
    a = np.array(a)
    mask = np.hstack((False, a[1:] > a[:-1])) & np.hstack((a[:-1] > a[1:], False)) & np.hstack((False, a[1:] > thresh))
    indices = np.arange(a.shape[0])
    return indices[mask], a[mask]

def getLocalMaxima(a):
    a = np.array(a)
    mask = np.hstack((False, a[1:] > a[:-1])) & np.hstack((a[:-1] > a[1:], False))
    indices = np.arange(a.shape[0])
    return indices[mask], a[mask]


def getGMMLastClusterCut(data):

    data = np.array(data).reshape((len(data), 1))

    # nComps = range(1, 3)
    # bics = []
    # aics = []
    # gmms = []
    #
    # for nComp in nComps:
    #     gmm = GaussianMixture(n_components=nComp)
    #     gmm.fit(data.reshape((data.shape[0], 1)))
    #     gmms.append(gmm)
    #     bics.append(gmm.bic(data))
    #     aics.append(gmm.aic(data))
    #
    # bestInd = np.argmin(bics)
    # bestGMM = gmms[bestInd]

    nComp = 5
    bestGMM = GaussianMixture(n_components=nComp)
    bestGMM.fit(data.reshape((data.shape[0], 1)))

    meansArgSort = np.argsort(bestGMM.means_[:, 0])
    means = bestGMM.means_[meansArgSort, 0]
    stds = np.sqrt(bestGMM.covariances_)[meansArgSort, 0].reshape((nComp,))

    hist, bins = np.histogram(data, np.arange(min(data), max(data), 0.1))
    hist = np.array(hist, dtype=float)

    indComps = map(lambda x: norm(loc=x[0], scale=x[1]), zip(means, stds))
    indDists = [x.pdf(bins) for x in indComps]
    indDists = [x / x.sum() for x in indDists]
    weights = bestGMM.weights_[meansArgSort]
    weightsDiagonalMatrix = np.diag(weights)
    weightedComps = np.dot(weightsDiagonalMatrix, np.vstack(tuple(indDists)))
    pdf = weightedComps.sum(axis=0)

    cdf = np.cumsum(pdf)

    # https://stackoverflow.com/questions/22579434/python-finding-the-intersection-point-of-two-gaussian-curves
    # intersection(s) between two Gaussians
    def intersectionsBetween2Gaussians(m1, m2, std1, std2):
        a = 1 / (2 * std1 ** 2) - 1 / (2 * std2 ** 2)
        b = m2 / (std2 ** 2) - m1 / (std1 ** 2)
        c = m1 ** 2 / (2 * std1 ** 2) - m2 ** 2 / (2 * std2 ** 2) - np.log(std2 / std1)
        return np.roots([a, b, c])

    # ind == left most cluster classified as spikes
    clusterCuts = []
    errors = []
    for ind in range(1, nComp):

        # pdf1 = weightedComps[:ind, :].sum(axis=0) / (weights[:ind].sum())
        # pdf2 = weightedComps[ind:, :].sum(axis=0) / (weights[ind:].sum())
        #
        # cdf1 = np.cumsum(pdf1)
        # cdf2 = np.cumsum(pdf2)
        #
        # mu1 = np.dot(pdf1, bins)
        # mu2 = np.dot(pdf2, bins)
        #
        # sig1 = np.dot(pdf1, np.power(bins - mu1, 2))
        # sig2 = np.dot(pdf2, np.power(bins - mu2, 2))

        mu1, mu2 = means[ind - 1], means[ind]
        sig1, sig2 = stds[ind - 1], stds[ind]

        if sig1 > 0 and sig2 > 0:
            # multiple intersections. choose the one between means
            clusterCutSuspects = intersectionsBetween2Gaussians(mu1, mu2, sig1, sig2)
            isBetweenMeans = map(lambda x: mu1 < x < mu2, clusterCutSuspects)
            if any(isBetweenMeans):
                clusterCutInd = np.where(isBetweenMeans)[0][0]
                clusterCut = clusterCutSuspects[clusterCutInd]
                closestBinInd = np.argmin(np.abs(bins - clusterCut))
                # discard cluster cuts which, when used for classification, result in less than 10% of identified maxima
                # being classified as spikes. This is done mainly to avoid large cluster cuts than can be caused by
                # far outliers.
                if cdf[closestBinInd] < 0.9:
                    error = mu2 - mu1
                    errors.append(error)
                    clusterCuts.append(clusterCut)

    if not clusterCuts:
        clusterCut = means[-1] - 3 * stds[1]
    else:
        clusterCut = clusterCuts[int(np.argmax(errors))]

    # plt.ion()
    # fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 5.6), sharex=True)
    #
    # for cc, er in zip(clusterCuts, errors):
    #     print("Cluster Cut={:.6f}, error={:.6f}".format(cc, er))
    # print("Means: {}".format(means))
    # print("Stds: {}".format(stds))
    # print("Mean Gaps: {}".format(np.diff(means)))
    #
    # ax[0].bar(bins[:-1], hist / hist.sum(), width=bins[1] - bins[0], color='b', label='data histogram')
    # ax[0].plot(bins[:-1], pdf[:-1], 'r', label='GMM fit')
    # ax[0].plot([clusterCut, clusterCut], ax[0].get_ylim(), 'r--', label='Cluster cut')
    # for mean in means:
    #     ax[0].plot([mean, mean], ax[0].get_ylim(), 'g--', label=None)
    # ax[0].plot([mean, mean], ax[0].get_ylim(), 'g--', label="Cluster means")
    #
    # # ax[0].set_xlabel('suspected spike  Heights')
    # ax[0].set_ylabel('Relative frequency')
    # ax[0].legend(loc='upper right')
    #
    # ax[1].bar(bins[:-1], np.cumsum(hist / hist.sum()), width=bins[1] - bins[0], color='b', label='data histogram')
    # ax[1].plot(bins[:-1], np.cumsum(pdf[:-1] / pdf[:-1].sum()), 'r', label='GMM fit')
    # ax[1].plot([clusterCut, clusterCut], ax[1].get_ylim(), 'r--', label='Cluster cut')
    # ax[1].set_xlabel('suspected spike  Heights')
    # ax[1].set_ylabel('CDF')
    # ax[1].legend(loc='upper right')
    #
    # fig.tight_layout()
    # raw_input()
    # plt.close(fig.number)

    return clusterCut


def filterButterworth(analogSignal, cutoff=100, transitionWidth=40, rippleDB=20):

    cutoff *= qu.Hz
    nyqFreq = analogSignal.sampling_rate / 2

    # b, a = iirfilter(order, cutoff / nyqFreq, btype='highpass', ftype='butter')

    transitionWidth = transitionWidth * qu.Hz

    N, beta = kaiserord(rippleDB, transitionWidth / nyqFreq)

    tapsLP = firwin(N, cutoff / nyqFreq, window=('kaiser', beta))

    temp = np.zeros((N,))
    temp[(N - 1) / 2] = 1
    tapsHP = temp - tapsLP

    delay = (N - 1) * 0.5 * analogSignal.sampling_period

    w, h = freqz(tapsHP, 1.0, int(nyqFreq))

    f = w * (nyqFreq) / (2 * PI)

    sigM = analogSignal.magnitude
    sigM = sigM.reshape((sigM.shape[0],))
    sigMFiltered = lfilter(tapsHP, 1.0, sigM)

    temp = AnalogSignal(
        signal=sigMFiltered,
        sampling_rate=analogSignal.sampling_rate,
        units=analogSignal.units,
        t_start=analogSignal.t_start - delay
    )
    filteredSignal = temp.reshape((temp.shape[0],))


    # from scipy.fftpack import fft
    # plt.ion()
    # fig1 = plt.figure()
    #
    # plt.plot(f, h, 'b-o')
    # plt.xlabel('Frequency(Hz)')
    # plt.ylabel('Gain')
    # plt.draw()
    #
    # fig2 = plt.figure()
    #
    #
    # plt.plot(analogSignal.times, analogSignal, 'r')
    # plt.plot(filteredSignal.times, filteredSignal, 'b')
    # plt.legend(['Original Signal', 'Filtered Signal'])
    # plt.xlabel('Time in ' + analogSignal.sampling_period.dimensionality.string)
    # plt.draw()
    #
    # fig3 = plt.figure()
    #
    #
    # freqResolution = 5 * qu.Hz
    # NFFT = np.ceil(2 * nyqFreq / freqResolution).simplified
    # NFFT = int(np.ceil(NFFT / 2)  * 2)
    # actualFreqResolution = 2 * nyqFreq / NFFT
    # frequencies = np.arange(int(NFFT/2)) * actualFreqResolution
    #
    #
    # outSig = filteredSignal.magnitude
    #
    #
    # origFFT = fft(sigM, NFFT)
    # origFFT[0] = 0
    # plt.plot(frequencies, np.abs(origFFT)[:int(NFFT/2)])
    #
    # filtFFT = fft(outSig, NFFT)
    # plt.plot(frequencies, np.abs(filtFFT)[:int(NFFT/2)])
    #
    # plt.legend(['FFT of Original Signal', 'FFT of Filtered Signal'])
    # plt.xlabel('Frequency(Hz)')
    # plt.draw()
    # raw_input()
    # plt.close(fig1.number)
    # plt.close(fig2.number)
    # plt.close(fig3.number)

    return filteredSignal


def getNoiseStd(analogSignal):
        sigMabs = np.abs(analogSignal.magnitude)

        # to remove regions of no signal, or signal artifically inserted by exclusion using " Intervals to Exclude (s)".
        sigMabsFiltered = sigMabs[sigMabs > 0.5]

        return np.median(sigMabsFiltered) / 0.6745 * analogSignal.units

        # return np.std(analogSignal.magnitude)


def thresholdAndDetect(analogSignal, minSpikeWidth, maxSpikeWidth):

    # include all maxima of noise as well potential spikes. Otherwise we might loose spikes as well.
    thres = getNoiseStd(analogSignal)
    print("Threshold={}".format(thres))

    maxSpikeWidthSamples = (maxSpikeWidth * qu.ms / analogSignal.sampling_period).simplified.magnitude

    minSpikeWidthSamples = (minSpikeWidth * qu.ms / analogSignal.sampling_period).simplified.magnitude

    sigM = analogSignal.magnitude
    sigM = sigM.reshape((sigM.shape[0],))
    thresholded = np.int8(sigM > thres)

    thresholdedDiff = np.diff(thresholded)

    onsets = np.where(thresholdedDiff > 0)[0] + 1

    offsets = np.where(thresholdedDiff < 0)[0] + 1

    spikeInds = []
    spikeHeights = []

    for onset, offset in zip(onsets, offsets):

        durAboveThesh = offset - onset
        if (durAboveThesh < maxSpikeWidthSamples) and (durAboveThesh >= minSpikeWidthSamples):
            spike = analogSignal[onset:offset + 1]
            spikeHeight = spike.max().magnitude
            spikeHeights.append(spikeHeight)
            spikeInd = (spike.argmax() + onset)
            spikeInds.append(spikeInd)

    spikeTimes = analogSignal.t_start + np.array(spikeInds) * analogSignal.sampling_period
    spikeHeights = np.array(spikeHeights) * analogSignal.units

    clusterCut = getGMMLastClusterCut(spikeHeights.magnitude)

    spikeTimes = spikeTimes[spikeHeights >= clusterCut]
    spikeHeights = spikeHeights[spikeHeights >= clusterCut]

    # fig4 = plt.figure()

    # plt.show(block=False)
    #
    # plt.plot(self.analogSignal.times, self.analogSignal, 'r')
    # plt.plot(self.filteredSignal.times, self.filteredSignal, 'b')
    # plt.plot(self.filteredSignal.times, np.ones_like(self.filteredSignal) * thres, 'g')
    # plt.legend(['Original Signal', 'Filtered Signal', 'Threshold'])
    # plt.plot(self.spikeTimes, self.spikePeaks, 'ro', ms=10)
    # plt.xlabel('Time in ' + self.analogSignal.sampling_period.dimensionality.string)
    # plt.draw()
    # raw_input("Press any key to continue...")

    return spikeTimes, spikeHeights


def detectSpikes(analogSignal, spikeMinWidth, spikeMaxWidth, *filterBW_args):

    filteredSignal = filterButterworth(analogSignal, *filterBW_args)
    signMed, filterSignalMedianRemoved = remove_median(filteredSignal)
    return thresholdAndDetect(filterSignalMedianRemoved, spikeMinWidth, spikeMaxWidth)


def remove_median(analogSignal):

    signMedian = np.median(analogSignal.magnitude) * analogSignal.units
    print("Signal Median={}".format(signMedian))
    return signMedian, analogSignal - signMedian


def removeSpikes(signal, template, enable_debugging):

    if enable_debugging:
        import ipdb

    signalMedian = signal - np.median(signal)
    templateMean = template - template.mean()

    corr = np.correlate(signalMedian.magnitude, templateMean.magnitude, mode='valid') / ((templateMean ** 2).sum())

    maxInds, maxs = getLocalMaxima(corr)

    if maxs.shape[0] == 0:

        return signal

    else:

        minInterval = int((3 * qu.ms * signal.sampling_rate).simplified)

        mask1 = (maxs > 0.5)

        maxs = maxs[mask1]
        maxInds = maxInds[mask1]

        if maxs.shape[0] == 0:
             return signal

        else:

            mask = [True]

            for ind in range(1, len(maxInds)):

                if (maxInds[ind] - maxInds[ind - 1]) < minInterval:

                    if maxs[ind] > maxs[ind - 1]:

                        mask[-1] = False
                        mask.append(True)

                    if maxs[ind] <= maxs[ind - 1]:

                        mask.append(False)

                else:
                    mask.append(True)

            mask = np.array(mask)

            maxs = maxs[mask]
            maxInds = maxInds[mask]

            if maxs.shape[0] == 0:
                return signal

            else:

                if enable_debugging:
                    fig, ax = plt.subplots()
                    ax.plot(corr, 'b')
                    ax.stem(maxInds, np.ones_like(maxInds), linefmt='r-o')
                    ax.set_xlabel('lags')
                    ax.set_ylabel('corr')

                    ipdb.set_trace()
                    plt.close(fig.number)

                toSub = np.zeros_like(signal)
                templateMin = template.min()
                for maxInd, maxVal in zip(maxInds, maxs):

                    toSub[maxInd: maxInd + template.shape[0]] += (template - templateMin)


                signal2Return = signalMedian - toSub + np.median(signal)

                if enable_debugging:
                    fig, ax = plt.subplots()
                    ax.plot(signal.times, signal, 'b')
                    ax.plot(signal2Return.times, signal2Return, 'r')
                    plt.draw()

                    ipdb.set_trace()
                    plt.close(fig.number)

                return signal2Return




