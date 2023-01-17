"""Console script for bnp_macs2."""
import typer
import dataclasses
from scipy.special import pdtrc
import numpy as np
import logging

from bionumpy.datatypes import Interval
from .control_pileup import get_control_pileup
from .fragment_pileup import get_fragment_pileup
# import matplotlib.pyplot as plt
import bionumpy as bnp

logging.basicConfig(level=logging.INFO)


@dataclasses.dataclass
class Params:
    fragment_length: int = 150
    p_value_cufoff: float = 0.001


# def get_control_pileup(reads, size, window_sizes, read_rate):
#     mid_points = (reads.start+reads.stop)//2
#     pileup = float(read_rate)
#     for window_size in window_sizes:
#         start = np.maximum(mid_points-window_size//2, 0)
#         stop = np.minimum(mid_points+window_size//2, size)
#         windows = dataclasses.replace(reads, start=start, stop=stop)
#         new_pileup = get_pileup(windows, size)
#         pileup = np.maximum(pileup, new_pileup/window_size)
#     return pileup


def logsf(count, mu):
    return np.log(pdtrc(count, mu))


def get_p_values(intervals, geometry, fragment_length, read_rate):
    fragment_pileup = get_fragment_pileup(intervals, fragment_length, chrom_size)
    control = fragment_length*get_control_pileup(intervals, [1000, 10000], read_rate, geometry)
    p_values = logsf(fragment_pileup, control)
    return p_values


@bnp.streamable()
def macs2(intervals, chrom_size, fragment_length, read_rate, p_value_cutoff, min_length, max_gap=30):
    if not len(intervals):
        return Interval.empty()
    p_values = get_p_values(intervals, chrom_size,
                            fragment_length, read_rate)
    peaks = call_peaks(p_values, p_value_cutoff, fragment_length)
    return Interval(intervals.chromosome[:len(peaks)], peaks.start, peaks.stop)


def main(filename: str, genome_file: str, fragment_length: int = 150, p_value_cutoff: float = 0.001, outfilename: str = None):
    genome = bnp.open(genome_file, buffer_type=bnp.io.files.ChromosomeSizeBuffer).read()
    genome_size = genome.size.sum()
    chrom_sizes = {str(name): size for name, size in zip(genome.name, genome.size)}
    intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read_chunks()
    multistream = bnp.MultiStream(chrom_sizes, intervals=intervals)
    n_reads = bnp.count_entries(filename)
    read_rate = n_reads/genome_size
    result = macs2(multistream.intervals, multistream.lengths, fragment_length,
                   read_rate, p_value_cutoff, fragment_length)
    if outfilename is not None:
        with bnp.open(outfilename, 'w') as f:
            f.write(result)
    return result


def test():
    #res = main("example_data/small_interval.bed", "example_data/small_genome.fa.fai")
    res = main('example_data/simulated_chip_seq.bed', 'example_data/simulated.chrom.sizes', fragment_length=200)
    for chunk in res:
        print(chunk)


def big():
    with bnp.open("/home/knut/Data/peaks.bed", "w") as outfile:
        for data in main("/home/knut/Data/subset.bed", "/home/knut/Data/hg38.chrom.sizes"):
            if len(data):
                outfile.write(data)


if __name__ == "__main__":

    typer.run(main)

def main():
    return


if __name__ == "__main__":
    typer.run(main)
