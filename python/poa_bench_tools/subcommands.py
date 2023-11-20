import argparse
import sys
import textwrap
from abc import ABCMeta, abstractmethod
from typing import Type


class Command(metaclass=ABCMeta):
    @classmethod
    @abstractmethod
    def register_arguments(cls, parser: argparse.ArgumentParser):
        pass

    @classmethod
    @abstractmethod
    def main(cls, args: argparse.Namespace):
        pass

    @classmethod
    def get_command_help(cls):
        subcommand_doc = cls.__doc__

        if not subcommand_doc:
            subcommand_doc = ""

        subcommand_doc = textwrap.dedent(subcommand_doc)
        first_help_line = subcommand_doc.strip().split('\n\n')[0].strip()

        return first_help_line, subcommand_doc


def run_command(command: Type[Command]):
    parser = argparse.ArgumentParser(description=command.__doc__)

    command.register_arguments(parser)

    args = parser.parse_args()
    sys.exit(int(command.main(args)))
