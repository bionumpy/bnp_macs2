from functools import reduce
from scipy.special import pdtrc as poisson_sf
import bionumpy as bnp
import numpy as np


def macs2(reads, params):
    fragment_pileup = reads.extended_to_size(params.fragment_length).get_pileup()
    local_averages = (reads.get_location('start').get_windows(window_size=w).get_pileup()/w
                      for w in params.window_sizes)
    global_average = params.n_reads/params.effective_genome_size
    control_pileup = reduce(np.maximum, local_averages, global_average) * params.fragment_length
    p_values = poisson_sf(fragment_pileup, control_pileup)
    peaks = bnp.GenomicIntervals.from_track(p_values < params.p_value_cutoff).merged(distance=params.max_gap)
    peaks = peaks[peaks.stop-peaks.start >= params.fragment_length]
    peak_signals = np.log10(p_values[peaks])*-1
    peaks, max_values, mean_values = bnp.compute((peaks, peak_signals.max(axis=-1), peak_signals.mean(axis=-1)))
    return bnp.datatypes.NarrowPeak(peaks.chromosome, peaks.start, peaks.stop,
                                    [f'peak_{i+1}' for i in range(len(peaks))],
                                    (max_values*10).astype(int),
                                    ['.']*len(peaks),
                                    mean_values,
                                    max_values,
                                    max_values,
                                    np.zeros_like(max_values, dtype=int))                                    
