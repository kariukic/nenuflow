#!/usr/bin/env python3

import os
import logging
from glob import glob
from pathlib import Path
import numpy as np
from casacore import tables as tab


import logging

logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.DEBUG)

from argparse import ArgumentParser

parser = ArgumentParser(description="Merge frequency subbands and split in time")

parser.add_argument(
    "-m",
    "--mslist",
    help="List of mses to be concatenated in frequency",
    nargs="+",
    required=True,
)

parser.add_argument(
    "-p",
    "--parset",
    help="DP3 parset",
    type=str,
    required=False,
)

parser.add_argument(
    "-n",
    "--nodes",
    help="List of nodes where the data is to be distributed",
    nargs="+",
    required=True,
)

parser.add_argument(
    "-o",
    "--msout",
    help="common name for the output MS files, each output MS will be labelled with a timechunk number added to this name",
    type=str,
    required=True,
)

parser.add_argument(
    "-d",
    "--datapath",
    help="The path where the data will be saved. Same for each node given",
    required=False,
    type=str,
)

parser.add_argument(
    "-c",
    "--datacolumn",
    help="The datacolumn to be used",
    required=False,
    type=str,
)

parser.add_argument(
    "-t",
    "--tolerance",
    default=0.25,
    help="If the last timechunk is longer than this fraction of the requested chuncklength (ntimes*ms_timestep), make it a separate MS. Otherwise add that chunk to the second last MS",
    required=False,
)

parser.add_argument(
    "-x",
    "--nmses_per_node",
    type=int,
    help="number of MS files to distribute to each node. The last node gets the remainder",
    required=False,
)

parser.add_argument(
    "-i", "--ntimes", type=int, help="Number of timesteps per output MS", required=False
)

parser.add_argument(
    "-q",
    "--nchunks",
    type=int,
    help="number of chunks to output instead of all possible ones based on ntims",
    required=False,
)

parser.add_argument(
    "-r",
    "--dry_run",
    help="Do a dry run first without running DP3",
    action="store_true",
)


def make_dir(dir_path) -> None:
    """
    Create a directory at the specified path if it does not already exist.

    Parameters:
    - dir_path (str): The path of the directory to be created.

    Returns:
    - None

    Example:
    >>> make_dir("/path/to/directory")
    """
    if not Path(dir_path).exists():
        try:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
            logging.debug(f"Made output directory: {dir_path}")
        except OSError as e:
            logging.error(f"Failed to create output directory: {dir_path}. Error: {e}")


def chunks(lst, n):
    """
    Yield successive n-sized chunks from a list.

    Parameters:
    - lst (list): The list from which to yield chunks.
    - n (int): The size of each chunk.

    Returns:
    - generator: A generator that yields successive n-sized chunks of the input list.

    Example:
    >>> lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    >>> for chunk in chunks(lst, 3):
    ...     print(chunk)
    [1, 2, 3]
    [4, 5, 6]
    [7, 8, 9]
    [10]
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def get_ms_duration(MS):
    """
    Calculate the total duration and timestep of a Measurement Set (MS) file.

    Parameters:
    - MS (str): The path to the MS file for which to calculate the duration and timestep.

    Returns:
    - tuple: A tuple containing the total duration (in seconds) and the timestep (in seconds) of the MS file.

    Example:
    >>> duration, timestep = get_ms_duration("/path/to/MS_file")
    """
    if os.path.isdir(MS):
        myMS = tab.table(MS)
    else:
        logging.error("Do not understand the format of MS", MS, "bailing out")
        return

    timerange = [
        np.amin(myMS.getcol("TIME_CENTROID")),
        np.amax(myMS.getcol("TIME_CENTROID")),
    ]
    timestep = myMS.getcell("INTERVAL", 0)

    logging.info(f"MS timerange: {timerange} s")
    logging.info(f"MS timestep: {timestep} s")

    duration = timerange[1] - timerange[0]  # / timestep
    logging.info(f"MS total duration: {duration/3600} hrs")

    return duration, timestep


def distribute(
    msout_name,
    starttimeslots,
    nodes,
    datapath,
    nmses_per_node,
    dry_run=False,
):
    """
    Distribute the output Measurement Sets (MS) across different nodes based on the specified parameters.

    Parameters:
    - msout_name (str): The base name for the output MS files.
    - starttimeslots (list): List of start time slots for the MS files.
    - nodes (list): List of node numbers to distribute the MS files.
    - datapath (str): The base data path where the MS files will be stored.
    - nmses_per_node (int): Number of MS files to be assigned to each node. If 0, the distribution will be based on the number of nodes.

    Returns:
    - list: A list of paths to the distributed MS files.

    Raises:
    - AssertionError: If the datapath does not start with '/data/' or if the number of starttimeslots does not match the number of generated MS file paths.

    Example:
    >>> distribute("output_MS", [0, 1, 2], [1, 2, 3], "/data/storage", 2)
    ['/net/node1/data/storage/output_MST001.MS', '/net/node1/data/storage/output_MST002.MS', '/net/node2/data/storage/output_MST003.MS', '/net/node2/data/storage/output_MST004.MS', '/net/node3/data/storage/output_MST005.MS', '/net/node3/data/storage/output_MST006.MS']
    """

    if msout_name[:-2] == "MS":
        msout = msout_name.strip("MS")
    if msout_name[:-2] == "ms":
        msout = msout_name.strip("ms")

    assert datapath.startswith("/data/")

    ms_nums = [t for t in range(len(starttimeslots))]
    msout_names = [msout_name + f"_T{t:03}.MS" for t in ms_nums]

    if nmses_per_node:
        msout_names = list(chunks(msout_names, nmses_per_node))
    else:
        msout_names = list(
            chunks(msout_names, int(np.ceil(len(starttimeslots) / len(nodes))))
        )

    logging.info(
        f"output MS nodes distribution: {list(zip([len(s) for s in msout_names], nodes))}"
    )

    all_msout_names = []
    for _n, (node, ms_sublist) in enumerate(zip(nodes, msout_names)):
        nodepath = f"/net/node{node}/{datapath}"
        if not dry_run:
            make_dir(nodepath)

        for ms in ms_sublist:
            all_msout_names.append(f"/net/node{node}/{datapath}/{ms}")

    assert len(starttimeslots) == len(all_msout_names)

    return all_msout_names


def splitMS(
    all_ms_files,
    datacolumn,
    ntimes,
    starttimeslots,
    all_msout_names,
    make_last_chunck_longer=False,
    parset=None,
    dry_run=False,
):
    for t, (tslot, msout) in enumerate(zip(starttimeslots, all_msout_names)):

        logging.debug(f"splitting from timeslot {tslot} to  {tslot+ntimes}")

        comm_prefix = f"DP3 {parset}" if parset else "DP3 steps=[]"

        if t == len(starttimeslots) - 1 and make_last_chunck_longer:
            comm = f"{comm_prefix} msin='{all_ms_files}' msin.datacolumn={datacolumn} msin.starttimeslot={tslot} msout.overwrite=true msout={msout}"
        else:
            comm = f"{comm_prefix} msin='{all_ms_files}' msin.datacolumn={datacolumn} msin.starttimeslot={tslot} msin.ntimes={ntimes} msout.overwrite=true msout={msout}"
        logging.info(f"DP3 command: {comm}")
        if not dry_run:
            os.system(comm)
        logging.info(f"Wrote timeslot {tslot} MS: {msout}")


def concatSubbands(
    all_ms_files,
    msout_name,
    datacolumn="DATA",
    ntimes=0,
    nchunks=0,
    nodes=None,
    datapath=None,
    nmses_per_node=1,
    tolerance=0.25,
    parset=None,
    dry_run=False,
):
    logging.info(f"Total MS files: {len(all_ms_files)}")
    if ntimes:
        duration, timestep = get_ms_duration(all_ms_files[0])

        starttimeslots = np.arange(0, duration // timestep, ntimes)

        if nchunks:
            assert nchunks < len(
                starttimeslots
            ), f"Number of time-chuncks requested {nchunks} is greater than total chuncks possible ({len(starttimeslots)}). Reduce `nchunks` or `ntimes`."
            starttimeslots = starttimeslots[:nchunks]

        remainder = (duration // timestep) % ntimes
        make_last_chunck_longer = False
        if remainder < tolerance * timestep * ntimes:
            starttimeslots = starttimeslots[:-1]
            make_last_chunck_longer = True
            logging.warning(f"The last MS will be of {remainder}s longer duration")
        else:
            logging.warning(
                f"The last MS will be of shorter duration: {remainder*timestep}s instead of {ntimes*timestep}s"
            )

    else:
        starttimeslots = [0]

    logging.debug(f"starttimeslot: {starttimeslots}")
    logging.info(f"total output MSes: {len(starttimeslots)}")

    if nodes:
        ntimes = int(ntimes)
        starttimeslots = [int(s) for s in starttimeslots]
        all_msout_names = distribute(
            msout_name, starttimeslots, nodes, datapath, nmses_per_node, dry_run=dry_run
        )
        splitMS(
            all_ms_files,
            datacolumn,
            ntimes,
            starttimeslots,
            all_msout_names,
            make_last_chunck_longer=make_last_chunck_longer,
            parset=parset,
            dry_run=dry_run,
        )


if __name__ == "__main__":
    args = parser.parse_args()
    all_ms_files = []
    for ms in args.mslist:
        if "*" in ms or "?" in ms:
            all_ms_files += glob(ms)
    else:
        all_ms_files = args.mslist

    try:
        assert all([os.path.isdir(ms) for ms in all_ms_files])
    except:
        all_ms_files = all_ms_files[0].strip().split(" ")
        print(all_ms_files[0])
        assert all([os.path.isdir(ms) for ms in all_ms_files])

    if args.parset:
        assert os.path.isfile(args.parset)

    concatSubbands(
        all_ms_files,
        args.msout,
        datacolumn=args.datacolumn,
        ntimes=args.ntimes,
        nchunks=args.nchunks,
        nodes=args.nodes,
        datapath=args.datapath,
        nmses_per_node=args.nmses_per_node,
        tolerance=args.tolerance,
        parset=args.parset,
        dry_run=args.dry_run,
    )
