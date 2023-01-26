"""Console script for bnp_macs2."""
import typer
import numpy as np
import logging

import bionumpy as bnp
from bionumpy.arithmetics.geometry import Geometry, StreamedGeometry, GenomicIntervals
from .macs2 import Macs2, Macs2Params
from .listener import Macs2Listner

logging.basicConfig(level=logging.INFO)


def main(filename: str,
         genome_file: str,
         fragment_length: int = 150,
         p_value_cutoff: float = 0.001,
         outprefix: str = None):
    genome = bnp.open(genome_file, buffer_type=bnp.io.files.ChromosomeSizeBuffer).read()
    chrom_sizes = {str(name): size for name, size in zip(genome.name, genome.size)}
    intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read()
    tag_size = np.median(intervals.stop-intervals.start)
    listner = Macs2Listner(lambda name: outprefix+name)
    params = Macs2Params(
        fragment_length=fragment_length,
        p_value_cutoff=p_value_cutoff,
        max_gap=int(tag_size),
        n_reads=bnp.count_entries(filename))

    m = Macs2(Geometry(chrom_sizes), params, listner)
    genomic_intervals = GenomicIntervals.from_intervals(intervals, chrom_sizes)
    return m.run(genomic_intervals)
    # genome_intervals = geometry.split_chromosomes(intervals)
    result = (geometry.apply(run, i).run(i) for name, i in groupby(intervals, 'chromosome'))
    return result


def run():
    typer.run(main)


if __name__ == "__main__":
    run()
