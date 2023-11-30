#!/usr/bin/env python3
"""
Sort a FASTA file by a Mash guide tree.

This script estimates pairwise genetic distances between sequences
in the FASTA file using mash, creates a NJ tree, and then
outputs a new FASTA file, where the sequences closest to each other
are output first.

This script indexes the FASTA with samtools faidx, which
should be available in your PATH.
"""

import argparse
import bz2
import gzip
import io
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam
import skbio
import numpy

from .subcommands import Command, run_command


def phylip_to_matrix(fname):
    with open(fname) as f:
        fiter = iter(f)
        num_entries = int(next(fiter))
        reference_names = []

        matrix = numpy.zeros((num_entries, num_entries))
        for i, line in enumerate(fiter):
            first_tab = line.find('\t')
            ref_name = line[:first_tab].strip()
            reference_names.append(ref_name)

            if first_tab >= 0:
                row = numpy.fromstring(line[first_tab+1:], sep='\t')
                matrix[i, 0:len(row)] = row

        index_upper = numpy.triu_indices(num_entries)
        matrix[index_upper] = matrix.T[index_upper]

    dmatrix = skbio.DistanceMatrix(matrix, ids=reference_names)

    return dmatrix


def create_tree_mash(fname: Path, k: int):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        subprocess.run(["mash", "sketch", "-i", "-k", str(k), fname, "-o", tmppath / "sketch"],
                       check=True)

        with open(tmppath / "dists.phylip", "wb") as ofile:
            subprocess.run(["mash", "triangle", tmppath / "sketch.msh"],
                           stdout=ofile, check=True)

        print("Creating tree with neighbor joining...", file=sys.stderr)
        dmatrix = phylip_to_matrix(tmppath / "dists.phylip")
        tree = skbio.tree.nj(
            dmatrix,
            result_constructor=lambda x: skbio.TreeNode.read(io.StringIO(x), convert_underscores=False)
        )

    return tree


class SortFasta(Command):
    @classmethod
    def register_arguments(cls, parser: argparse.ArgumentParser):
        parser.add_argument(
            '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help="Output filename for the sorted FASTA, defaults to stdout."
        )

        parser.add_argument(
            '-t', '--tree', type=Path, default=None, required=False,
            help="Path to a guide tree in newick format. If not given, will use Mash to create one."
        )

        parser.add_argument(
            '-k', '--kmer_size', type=int, default=15,
            help="Mash k-mer size for sketching"
        )

        parser.add_argument(
            '-O', '--tree-output', type=argparse.FileType('w'), default=None,
            help="When creating a guide-tree using mash, save the tree in Newick format to the given file."
        )

        parser.add_argument('fasta', type=Path, metavar='FASTA', nargs='+',
                            help="Path to FASTA file to sort.")

    @classmethod
    def main(cls, args: argparse.Namespace):
        with tempfile.TemporaryDirectory() as tmpdir:
            print("Combining and indexing input FASTA files...", file=sys.stderr)
            tmppath = Path(tmpdir)

            with open(tmppath / "all_seq.fna", "wb") as ofile:
                for fasta in args.fasta:
                    open_func = {
                        '.gz': gzip.open,
                        '.bz2': bz2.open,
                    }.get(fasta.suffix, open)

                    with open_func(fasta, "rb") as ifile:
                        ofile.write(ifile.read())

            subprocess.run(["bgzip", tmppath / "all_seq.fna"], check=True)
            subprocess.run(["samtools", "faidx", tmppath / "all_seq.fna.gz"])

            if args.tree:
                with open(args.tree) as ifile:
                    tree = skbio.TreeNode.read(ifile)
            else:
                tree = create_tree_mash(tmppath / "all_seq.fna.gz", args.kmer_size)

                if args.tree_output is not None:
                    tree.write(args.tree_output, "newick")

            dmatrix = tree.tip_tip_distances()

            print("Writing sorted FASTA entries...", file=sys.stderr)
            reader = pysam.FastaFile(tmppath / "all_seq.fna.gz")
            for seq_id in dmatrix.ids:
                try:
                    seq = reader.fetch(seq_id)
                except ValueError:
                    print("Could not read sequence", seq_id, file=sys.stderr)
                    continue

                skbio.io.write(
                    skbio.DNA(seq, metadata={"id": seq_id}),
                    "fasta",
                    into=args.output
                )

            print("Done.", file=sys.stderr)


SortFasta.__doc__ = __doc__


if __name__ == '__main__':
    run_command(SortFasta)
