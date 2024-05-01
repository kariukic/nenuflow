#! /usr/bin/env python3
from argparse import ArgumentParser, FileType
import numpy as np
import sys
import os

parser = ArgumentParser("get source names from the modeltool attenuate log file")
parser.add_argument(
    "-i",
    "--model_attenuate_log",
    help="Log file produced by model attenuate",
    dest="model_attenuate_log",
    required=True,
    type=FileType('r'),
)

parser.add_argument(
    "-o",
    "--output_txt",
    help="output txt file",
    dest="output_txt",
    default="apparent_sources_left.txt",
    type=FileType('w'),
)

parser.add_argument (
    "-e",
    "--excude",
    nargs='+',
    help="Exclude this source in outputs",
    dest="exclude",
    default=[],
)

#we want to extract the source names from such a part in the modeltool attenuate log file
"""
Sky model info after model editing:
- CasA: 24 cmpts totaling 148.3 Jy
- CygA: 8 cmpts totaling 45.2 Jy
- Main: 758 cmpts totaling 1347.6 Jy
- VirA: 28 cmpts totaling 136.4 Jy
Saving apparent sky model to apparent_sky.catalog ...
"""


def main(argv):
    args = parser.parse_args(argv)
    log_file = args.model_attenuate_log

    lines = log_file.readlines()
    for l, ll in enumerate(lines):
        if ll.startswith("Sky model info after model editing:"):
            start = l
        if ll.startswith("Saving apparent sky model to"):
            stop = l
    sources = [s.split('- ')[1].split(':')[0] for s in lines[start+1:stop]]
    sources = [[s] for s in sources if s not in args.exclude]

    # with open(args.output_txt, mode='w') as a:
    args.output_txt.write(str(sources).replace(" ", ""))



if __name__ == "__main__":
    main(sys.argv[1:])
