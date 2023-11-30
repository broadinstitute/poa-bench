#!/usr/bin/env python3
import gzip
import argparse
import itertools
from pathlib import Path

import numpy.random

from .create_mutated_sequences import create_mutated_sequences
from .subcommands import Command, run_command


NUM_GRAPH_HAP = [3, 10, 25]
ERROR_RATES = [0.01, 0.03, 0.05]
LENGTHS = [500, 1000, 1500]


class MakeSyntheticData(Command):
    """
    Create multiple synthetic benchmark sets by randomly mutating
    a set of base sequences.
    """

    @classmethod
    def register_arguments(cls, parser: argparse.ArgumentParser):
        parser.add_argument(
            'base', type=Path,
            help="FASTA file with a single sequence used as base for creating mutated sequences."
        )

        parser.add_argument(
            '-H', '--num-haplotypes', type=int, nargs='+',
            help="Number of base haplotypes to generate."
        )

        parser.add_argument(
            '-E', '--error-rates', type=float, nargs='+',
            help="Error rates to use."
        )

        parser.add_argument(
            '-L', '--truncate-lengths', type=int, nargs='+',
            help="Truncate the input sequence."
        )

        parser.add_argument(
            '-N', '--num-sequences', type=int, default=100,
            help="Number of randomly mutated sequences to generate"
        )

        parser.add_argument(
            '-o', '--output-dir', type=Path,
            help="Output directory for each dataset"
        )

        parser.add_argument(
            '-s', '--seed', type=int, default=None, required=False,
            help="Random generator seed"
        )

    @classmethod
    def main(cls, args: argparse.Namespace):
        num_graph_hap = args.num_haplotypes if args.num_haplotypes else NUM_GRAPH_HAP
        error_rates = args.error_rates if args.error_rates else ERROR_RATES
        lengths = args.truncate_lengths if args.truncate_lengths else [None]

        rng = numpy.random.default_rng(args.seed)

        for num_hap, error_rate, length in itertools.product(num_graph_hap, error_rates, lengths):
            if length:
                output_dir = args.output_dir / f"hap{num_hap}_err{error_rate:g}_len{length}"
            else:
                output_dir = args.output_dir / f"hap{num_hap}_err{error_rate:g}"

            output_dir.mkdir(parents=True, exist_ok=True)

            # First, generate random haplotypes to be used for the graph
            graph_seq = output_dir / f"base_haplotypes.fna.gz"
            with gzip.open(graph_seq, "wt") as o:
                create_mutated_sequences(args.base, o, error_rate, num_hap, length, rng)

            # Then, using the just generated graph haplotypes as base, generate more random sequences to be used
            # for alignment benchmarking
            align_seq = output_dir / f"align_seq.fna.gz"
            with gzip.open(align_seq, "wt") as o:
                create_mutated_sequences(graph_seq, o, error_rate, args.num_sequences, length, rng)

            with open(output_dir / "meta.toml", "wt") as o:
                print(f"error_rate = {error_rate:g}", file=o)
                print(f"length = {length}", file=o)
                print(file=o)

                print("[graph_set]", file=o)
                print(f'fname = "{graph_seq.name}"', file=o)
                print(f"num_hap = {num_hap}", file=o)
                print(file=o)

                print("[align_set]", file=o)
                print(f'fname = "{align_seq.name}"', file=o)


if __name__ == '__main__':
    run_command(MakeSyntheticData)
