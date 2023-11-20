#!/usr/bin/env python3
"""
Main entry point for the poa_bench_tools command
"""

import argparse
from importlib.metadata import entry_points


def main():
    parser = argparse.ArgumentParser(
        description="Helper script to create and manage data sets to be used for benchmarking POA tools."
    )

    subparsers = parser.add_subparsers(help="Description", metavar='SUBCOMMAND')

    for ep in entry_points(group='poa_bench_tools.subcommands'):
        command = ep.load()
        first_help_line, description = command.get_command_help()

        subparser = subparsers.add_parser(ep.name, help=first_help_line, description=description)
        command.register_arguments(subparser)
        subparser.set_defaults(subcommand_func=command.main)

    args = parser.parse_args()
    args.subcommand_func(args)

    return 0
