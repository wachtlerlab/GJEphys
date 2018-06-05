from scipy.optimize import least_squares
import quantities as qu
import numpy as np

def estimateGJStimPars(stim):

    def modelGJStim(x, stim):

        [A, f, phi, B] = x

        stimTimes = stim.times
        stimTimes.unit = qu.s

        stimTimesS = stimTimes.magnitude

        modelStim = B + A * np.sin(2 * np.pi * f * stimTimesS + phi)

        return stim.magnitude - modelStim

    optRes = least_squares(modelGJStim, [10, 265, 0, 0], args=(stim,), bounds=([5, -np.inf, 0, -np.inf],
                                                                            [50, np.inf, 2 * np.pi, np.inf]))

    return optRes.x[:3]




