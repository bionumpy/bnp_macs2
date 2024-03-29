"""Console script for bnp_macs2."""
import typer
import numpy as np
import logging

import bionumpy as bnp
from bionumpy.genomic_data import Genome
from .macs2 import Macs2, Macs2Params
from .listener import Macs2Listner, StreamListner

logging.basicConfig(level=logging.INFO)

stream = False


def main(filename: str,
         genome_file: str,
         fragment_length: int = 150,
         p_value_cutoff: float = 0.001,
         outprefix: str = None):

    genome = Genome.from_file(genome_file)
    # bnp.open(genome_file, buffer_type=bnp.io.files.ChromosomeSizeBuffer).read()
    # chrom_sizes = {str(name): size for name, size in zip(genome.name, genome.size)}
    tmp = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read_chunk()
    tag_size = np.median(tmp.stop-tmp.start)
    if stream:
        intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read_chunks()
    else:
        intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read()
    # intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read_chunks()
    listner = Macs2Listner(lambda name: outprefix+name)
    listner = StreamListner(lambda name: outprefix+name)
    params = Macs2Params(
        fragment_length=fragment_length,
        p_value_cutoff=p_value_cutoff,
        max_gap=int(tag_size),
        n_reads=bnp.count_entries(filename),
        effective_genome_size = genome.size)

    m = Macs2(params, listner)
    intervals = genome.read_intervals(filename, stream=False, stranded=True)
    # genomic_intervals = genome.open(filename, ...) -> Ineterf
    #if stream:
    #    genomic_intervals = GenomicIntervals.from_interval_stream(intervals, chrom_sizes)
    #else:
    #     genomic_intervals = GenomicIntervals.from_intervals(intervals, chrom_sizes)
    return m.run(intervals)


def run():
    typer.run(main)


if __name__ == "__main__":
    run()
