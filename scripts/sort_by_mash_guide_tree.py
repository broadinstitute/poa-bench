#!/usr/bin/env python3
"""
Sort a FASTA file by a Mash guide tree.

This script estimates pairwise genetic distances between sequences
in the FASTA file using mash, creates a NJ tree, and then
outputs a new FASTA file, where the sequences closest to each other
are output first.

Requires that the FASTA file is indexed using `samtools faidx`.
"""

import argparse
import io
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam
import skbio
import numpy


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


def create_tree_mash(fname: Path):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        subprocess.run(["mash", "sketch", "-i", "-k", "11", fname, "-o", tmppath / "sketch"],
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


def main():
    parser = argparse.ArgumentParser(
        description='Sort FASTA by genetic distance using a guide tree'
    )

    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help="Output filename for the sorted FASTA, defaults to stdout."
    )

    parser.add_argument(
        '-t', '--tree', type=Path, default=None, required=False,
        help="Path to a guide tree in newick format. If not given, will use Mash to create one."
    )

    parser.add_argument('fasta', type=Path, metavar='FASTA',
                        help="Path to FASTA file to sort.")

    args = parser.parse_args()

    if args.tree:
        with open(args.tree) as ifile:
            tree = skbio.TreeNode.read(ifile)
    else:
        tree = create_tree_mash(args.fasta)

    dmatrix = tree.tip_tip_distances()

    print("Writing sorted FASTA entries...", file=sys.stderr)
    reader = pysam.FastaFile(args.fasta)
    for seq_id in dmatrix.ids:
        seq = reader.fetch(seq_id)
        skbio.io.write(
            skbio.DNA(seq, metadata={"id": seq_id}),
            "fasta",
            into=args.output
        )

    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()







