#!/usr/bin/env python3
"""
Analyze a multiple sequence alignment and oytput
various statistics (e.g., ANI, average length, ...)
"""

import argparse
import gzip
import itertools
import multiprocessing
import sys
from pathlib import Path

import skbio
import numpy

from .subcommands import Command, run_command


def compute_ani(seqs: tuple[numpy.array, numpy.array]):
    seq1, seq2 = seqs
    ani = numpy.char.equal(seq1, seq2).sum() / seq1.size

    return ani


class MSAStats(Command):
    @classmethod
    def register_arguments(cls, parser: argparse.ArgumentParser):
        parser.add_argument(
            'fasta', type=Path,
            help="Input MSA in FASTA format"
        )

        parser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help="Output file. Defaults to stdout."
        )

        parser.add_argument(
            '-j', '--parallel', type=int, default=None,
            help="Number of parallel processes to start"
        )

    @classmethod
    def main(cls, args: argparse.Namespace):
        open_func = {
            '.gz': gzip.open,
        }.get(args.fasta.suffix, open)

        seqs: list[numpy.array] = []
        with open_func(args.fasta, "rt") as ifile:
            for r in skbio.io.read(ifile, "fasta"):
                seqs.append(numpy.array([str(r).encode('ascii')], dtype=numpy.bytes_).view('S1').reshape((len(r),)))

        seq_lengths = numpy.array([
            (seq != b'-').sum() for seq in seqs
        ])

        print("num_seqs", len(seqs), sep='\t', file=args.output)
        print("seq_length_mean", seq_lengths.mean(), sep='\t', file=args.output)
        print("seq_length_std", seq_lengths.std(), sep='\t', file=args.output)

        # Compute pairwise ANI
        anis = []
        with multiprocessing.Pool(args.parallel) as p:
            anis = list(p.imap_unordered(
                compute_ani,
                itertools.product(seqs, repeat=2),
                chunksize=64,
            ))

        print("ani", numpy.array(anis).mean(), sep='\t', file=args.output)


MSAStats.__doc__ = __doc__


if __name__ == '__main__':
    run_command(MSAStats)
