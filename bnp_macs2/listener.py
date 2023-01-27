import bionumpy as bnp
from bionumpy.arithmetics.geometry import GenomicTrack
import matplotlib.pyplot as plt
import numpy as np


class Listner:
    def __init__(self, filename_template, names=[]):
        self._filename_template = filename_template
        self._names = names


class Macs2Listner(Listner):
    def control_lambda(self, track: GenomicTrack):
        plt.plot(track._global_track.to_array(), label='control_lambda'); plt.legend(); plt.show()
        bnp.open(self._filename_template('control_lambda.bdg'), 'w').write(track.to_bedgraph())

    def treat_pileup(self, track: GenomicTrack):
        plt.plot(track._global_track.to_array(), label='treat_pileup')
        bnp.open(self._filename_template('treat_pileup.bdg'), 'w').write(track.to_bedgraph())

    def p_scores(self, track):
        plt.plot(track._global_track.to_array(), label='p_scores'); plt.legend(); plt.show()

    def peaks(self, peaks: bnp.datatypes.NarrowPeak):
        peaks.signal_value[peaks.signal_value == np.inf] = 10000
        peaks.score[peaks.score == np.inf] = 1000
        peaks.score = peaks.score.astype(int)
        peaks.p_value[peaks.p_value == np.inf] = 10000
        peaks.q_value[peaks.q_value == np.inf] = 10000
        print(peaks)
        bnp.open(self._filename_template('peaks.narrowPeak'), 'w').write(peaks)

class register:
    def __init__(self, name):
        self._name = name

    def __call__(self, func):
        def new_func(obj, *args, **kwargs):
            result = func(obj, *args, **kwargs)
            if hasattr(obj, '_listner') and obj._listner is not None:
                getattr(obj._listner, self._name)(result)
            return result
        return new_func

