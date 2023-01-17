import numpy as np


def call_peaks(log_p_values, p_value_cufoff, min_length, max_gap, geometry):
    peaks = log_p_values < np.log(p_value_cufoff)
    peaks = peaks.to_intervals()
    peaks = geometry.merge_intervals(peaks, distance=max_gap)
    peaks = peaks[(peaks.stop-peaks.start) >= min_length]
    return peaks

