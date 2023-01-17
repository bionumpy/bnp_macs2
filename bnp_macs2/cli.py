"""Console script for bnp_macs2."""
import typer
import dataclasses
from scipy.special import pdtrc
import numpy as np
import logging

from bionumpy.datatypes import Interval
from bionumpy.arithmetics.geometry import Geometry
from .control_pileup import get_control_pileup
from .fragment_pileup import get_fragment_pileup
from .call_peaks import call_peaks
# import matplotlib.pyplot as plt
import bionumpy as bnp

logging.basicConfig(level=logging.INFO)


@dataclasses.dataclass
class Macs2Params:
    fragment_length: int = 150
    n_reads: int = None
    p_value_cutoff: float = 0.001
    max_gap: int = 30

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
    fragment_pileup = get_fragment_pileup(intervals, fragment_length, geometry)
    control = fragment_length*get_control_pileup(intervals, [1000, 10000], read_rate, geometry)
    p_values = logsf(fragment_pileup, control)
    return p_values


# @bnp.streamable()
def macs2(intervals: Interval, geometry: Geometry, params: Macs2Params):
    if not len(intervals):
        return Interval.empty()
    p_values = get_p_values(intervals, geometry,
                            params.fragment_length,
                            params.n_reads/geometry.size())
    return call_peaks(p_values, params.p_value_cutoff, params.fragment_length, params.max_gap, geometry)


def main(filename: str, genome_file: str, fragment_length: int = 150, p_value_cutoff: float = 0.001, outfilename: str = None):
    genome = bnp.open(genome_file, buffer_type=bnp.io.files.ChromosomeSizeBuffer).read()
    chrom_sizes = {str(name): size for name, size in zip(genome.name, genome.size)}
    intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read()
    tag_size = np.median(intervals.stop-intervals.start)
    params = Macs2Params(
        fragment_length=fragment_length,
        p_value_cutoff=p_value_cutoff,
        max_gap=int(tag_size),
        n_reads=bnp.count_entries(filename))
    result = macs2(intervals, Geometry(chrom_sizes), params)
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


def run():
    typer.run(main)


if __name__ == "__main__":
    typer.run(main)
