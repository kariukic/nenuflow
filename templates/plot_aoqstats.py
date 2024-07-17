import aoquality as ao
import numpy as np
import matplotlib.pyplot as plt

from argparse import ArgumentParser, FileType

parser = ArgumentParser(description="Plot AOquality stats")

parser.add_argument(
    "-q",
    "--qs_file",
    help="Aoquality statistics file",
    required=True,
)

parser.add_argument(
    "-s",
    "--stat",
    help="statistic to plot",
    default="Std",
    required=False,
    type=str,
)

parser.add_argument(
    "-o",
    "--output_file",
    help="output txt file",
    default="stats.png",
)

if __name__ == "__main__":
    args = parser.parse_args()
    aoq = ao.AOQualityTimeStat(args.qs_file)
    stat_data = aoq.get_stat(args.stat)
    time = aoq.time
    output_data = np.zeros((len(time), 3), dtype="complex")
    output_data[:, 0] = aoq.time
    output_data[:, 1] = stat_data[:, 0]
    output_data[:, 2] = stat_data[:, 3]
    data = output_data.real

    fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
    fig.set_size_inches((12, 6))
    axes[0].plot((data[:, 0] - data[:, 0][0]) / 3600, data[:, 1])
    axes[1].plot((data[:, 0] - data[:, 0][0]) / 3600, data[:, 2])
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"Std (real part of visibilities)")
    axes[0].set_xlabel("Time (hour)")
    axes[1].set_xlabel("Time (hour)")
    axes[0].set_title("XX")
    axes[1].set_title("YY")
    fig.tight_layout()
    fig.savefig(args.output_file, dpi=100)
