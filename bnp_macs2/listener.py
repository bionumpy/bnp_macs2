import bionumpy as bnp
from bionumpy.genomic_data import GenomicArray
import matplotlib.pyplot as plt
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Listner:
    pass


class FileListner(Listner):
    def __init__(self, filename_template, names=[]):
        self._filename_template = filename_template
        self._names = names


class Macs2Listner(FileListner):
    def control_lambda(self, track: GenomicArray):
        plt.plot(track._global_track.to_array(), label='control_lambda'); plt.legend(); plt.show()
        bnp.open(self._filename_template('control_lambda.bdg'), 'w').write(track.to_bedgraph())

    def treat_pileup(self, track: GenomicArray):
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
        bnp.open(self._filename_template('peaks.narrowPeak'), 'w').write(peaks)


class StreamListner(FileListner):
    def peaks(self, peaks: bnp.datatypes.NarrowPeak):
        peaks.signal_value[peaks.signal_value == np.inf] = 10000
        peaks.score[peaks.score == np.inf] = 1000
        peaks.score = peaks.score.astype(int)
        peaks.p_value[peaks.p_value == np.inf] = 10000
        peaks.q_value[peaks.q_value == np.inf] = 10000
        bnp.open(self._filename_template('peaks.narrowPeak'), 'w').write(peaks)


class register:
    def __init__(self, name):
        self._name = name

    def __call__(self, func):
        def new_func(obj, *args, **kwargs):
            result = func(obj, *args, **kwargs)
            logger.info(f'Registering {self._name}')
            if hasattr(obj, '_listner') and obj._listner is not None and hasattr(obj._listner, self._name):
                getattr(obj._listner, self._name)(result)
            return result
        return new_func


class DebugListnerStream(Listner):

    def treat_pileup(self, fragment_pileup):
        print('fragment pilup', fragment_pileup)
        # fragment_pileup.add_callback(print)
