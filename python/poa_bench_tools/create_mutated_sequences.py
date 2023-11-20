#!/usr/bin/env python3
"""
Generate synthetically mutated sequences at a given error rate
using a set of base haplotypes.
"""

import argparse
import bz2
import gzip
import math
import sys
from pathlib import Path

import numpy.random
import skbio
import numpy as np

from .subcommands import Command, run_command


events = ["S", "D", "I"]
p_event = [0.75, 0.125, 0.125]  # SNP, deletion, insertion
nucleotides = "ACGT"
nucleotides_set = set(nucleotides)


def open_compressed(fname: Path, *args, **kwargs):
    func = {
        '.gz': gzip.open,
        '.bz2': bz2.open,
    }.get(fname.suffix, open)

    return func(fname, *args, **kwargs)


def random_seq(length, rng=None) -> str:
    if rng is None:
        rng = np.random.default_rng()

    return "".join(str(v) for v in rng.choice(list(nucleotides), size=length))


def mutate_seq(seq, error_rate, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    num_errors = int(round(len(seq) * error_rate))
    if num_errors == 0:
        return seq

    mutations_pos = rng.choice(np.arange(len(seq)), num_errors, replace=False)
    mutations_pos.sort()

    components = []
    if mutations_pos[0] > 0:
        components.append(seq[:mutations_pos[0]])

    next_index = 0
    for i, pos in enumerate(mutations_pos):
        event = rng.choice(events, p=p_event)

        if pos < next_index:
            continue

        match event:
            case "S":
                remaining = list(nucleotides_set - {seq[pos]})
                sub = str(rng.choice(remaining))
                components.append(sub)
                next_index = pos + 1
            case "I":
                length = rng.integers(0, 3)
                insertion = random_seq(length, rng)
                components.append(seq[pos])
                components.append(insertion)
                next_index = pos + 1
            case "D":
                length = rng.integers(0, 3)
                next_index = pos + length

        end_pos = mutations_pos[i+1] if i + 1 < len(mutations_pos) else None
        components.append(seq[next_index:end_pos])

    return "".join(components)


def create_mutated_sequences(base_haplotypes: Path, output_file, error_rate: float, num: int,
                             truncate: int = None, rng: numpy.random.Generator = None):
    if rng is None:
        rng = numpy.random.default_rng()

    num_haplotypes = 0
    with open_compressed(base_haplotypes, "rt") as f:
        for _ in skbio.io.read(f, "fasta"):
            num_haplotypes += 1

    num_seq_per_hap = int(math.ceil(num / num_haplotypes))

    with open_compressed(base_haplotypes, "rt") as f:
        total_generated = 0
        for r in skbio.io.read(f, "fasta"):
            print("Mutating", r.metadata['id'], file=sys.stderr)

            for i in range(num_seq_per_hap):
                if total_generated >= num:
                    break

                sequence = str(r)

                if truncate is not None and truncate > 0:
                    sequence = sequence[:truncate]

                mutated = mutate_seq(sequence, error_rate, rng)
                seq = skbio.DNA(mutated, metadata={
                    "id": r.metadata['id'] + f"_mut{i}",
                    "description": r.metadata['description']
                })

                skbio.io.write(seq, "fasta", into=output_file)

                total_generated += 1


class MutateSeq(Command):
    @classmethod
    def register_arguments(cls, parser: argparse.ArgumentParser):
        parser.add_argument(
            'base_haplotypes', type=Path, metavar='BASE_HAP',
            help="FASTA file with base haplotypes to mutate"
        )

        parser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help="Output file to write mutated sequences to. Default: stdout."
        )

        parser.add_argument(
            '-e', '--error', type=float, default=0.01,
            help="Error rate. Default: %(default)g."
        )

        parser.add_argument(
            '-n', '--num-to-generate', type=int, default=1000,
            help="Number of mutated sequences to generate per base haplotype."
        )

        parser.add_argument(
            '-t', '--truncate', type=int, default=None, required=False,
            help="Truncate length of sequences to the given number. Optional."
        )

        parser.add_argument(
            '-s', '--seed', type=int, default=None, required=False,
            help="Random generator seed. Optional."
        )

    @classmethod
    def main(cls, args: argparse.Namespace):
        rng = np.random.default_rng(args.seed)

        create_mutated_sequences(args.base_haplotypes, args.output, args.error_rate, args.num_to_generate,
                                 args.truncate, rng)


MutateSeq.__doc__ = __doc__


if __name__ == '__main__':
    run_command(MutateSeq)
