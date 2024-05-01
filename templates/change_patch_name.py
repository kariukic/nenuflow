#! /usr/bin/env python3
from argparse import ArgumentParser, FileType
import numpy as np
import sys
import os

parser = ArgumentParser("Change the name of a patch")
parser.add_argument(
    "-i",
    "--model",
    help="Log file produced by model attenuate",
    dest="model",
    required=True,
)

parser.add_argument(
    "-x",
    "--patch_name_in",
    dest="patch_name_in",
    type=str,
    help="patch name in",
)

parser.add_argument (
    "-y",
    "--patch_name_out",
    dest="patch_name_out",
    type=str,
    help="patch name out",
)


def main(argv):
    args = parser.parse_args(argv)

    with open(args.model, mode = 'r') as mi:
        lines = [l.replace(args.patch_name_in, args.patch_name_out) for l in mi.readlines()]

    with open(args.model, mode = 'w') as mo:
        mo.writelines(lines)

open

if __name__ == "__main__":
    main(sys.argv[1:])