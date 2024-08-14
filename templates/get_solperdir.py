#!/usr/bin/env python3
# Author: S. Munshi

import numpy as np
import lsmtool
from argparse import ArgumentParser
import sys


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


parser = ArgumentParser(description="Return solutions per direction")
parser.add_argument("skymodel", help="Path to clustered apparent skymodel", type=str)
parser.add_argument(
    "--solint", help="Solint parameter in DP3", type=int, required=False, default=48
)
parser.add_argument(
    "--max_flux",
    help="Maximum flux in a single solution",
    type=float,
    required=False,
    default=40.0,
)
parser.add_argument(
    "--min_flux",
    help="Minimum flux in a single solution",
    type=float,
    required=False,
    default=15.0,
)
parser.add_argument(
    "--dist",
    help="Include angular distance of patch in solint calculation",
    type=str2bool,
    required=False,
    default=False,
)
args = parser.parse_args(sys.argv[1:])

# # Without spectral index based frequency correction
# def get_patch_fluxes(skymodel):
#     sm = lsmtool.load(skymodel)
#     patches = sm.getPatchNames()
#     patch_fluxes = sm.getColValues('I', aggregate='sum')
#     main = sm.getPatchPositions('Main')
#     separations = sm.getDistance(main['Main'][0].value,main['Main'][1].value,'Main')
#     return list(patches), list(patch_fluxes), list(separations)


def get_flux(flux, alpha, logSI, nterms, ref_freq, target_freq):
    if nterms > 1:
        for i in range(nterms):
            if logSI == "true":
                flux *= 10.0 ** (
                    alpha[i] * (np.log10(target_freq / ref_freq)) ** (i + 1)
                )
            else:
                flux += alpha[i] * ((target_freq / ref_freq) - 1.0) ** (i + 1)
    else:
        if logSI == "true":
            flux *= 10.0 ** (alpha[0] * np.log10(target_freq / ref_freq))
        else:
            flux += alpha[0] * ((target_freq / ref_freq) - 1.0)
    return flux


def get_patch_fluxes(skymodel, targetFreq=66.5e6):
    sm_full = lsmtool.load(skymodel)
    patches = sm_full.getPatchNames()
    patch_fluxes = []
    for i in range(len(patches)):
        sm = lsmtool.load(skymodel)
        sm.select("Patch == %s" % (patches[i]))
        if "ReferenceFrequency" in sm.getColNames():
            refFreq = sm.getColValues("ReferenceFrequency")
        else:
            refFreq = np.array([sm.table.meta["ReferenceFrequency"]] * len(sm))
        alphas = sm.getColValues("SpectralIndex")
        nterms = alphas.shape[1]
        fluxes = sm.getColValues("I")
        if "LogarithmicSI" in sm.getColNames():
            logSI = sm.getColValues("LogarithmicSI")
        else:
            logSI = np.array(["true"] * len(sm))
        for i in range(len(sm)):
            fluxes[i] = get_flux(
                fluxes[i], alphas[i], logSI[i], nterms, refFreq[i], targetFreq
            )
        patch_fluxes.append(sum(fluxes))
    main = sm.getPatchPositions("Main")
    separations = sm.getDistance(main["Main"][0].value, main["Main"][1].value, "Main")
    return list(patches), list(patch_fluxes), list(separations)


def get_divisors(n):
    divisors = []
    for i in range(1, n + 1):
        if n % i == 0:
            divisors.append(i)
    return divisors


def find_closest_element(lst, target):
    return min(lst, key=lambda x: abs(x - target))


def get_solints(skymodel, solint=48, max_flux=40, min_flux=15, dist=False):
    patches, patch_fluxes, separations = get_patch_fluxes(skymodel)
    main_idx = patches.index("Main")
    patch_fluxes.pop(main_idx)
    flux_norm = np.array(patch_fluxes) / np.maximum(
        np.minimum(max_flux, min(patch_fluxes)), min_flux
    )
    if dist:
        separations.pop(main_idx)
        vel_norm = np.array(separations) / min(separations)
        weights_norm = np.sqrt(vel_norm * flux_norm)
    else:
        weights_norm = flux_norm
    solints = get_divisors(solint)
    weights_solint = [find_closest_element(solints, i) for i in weights_norm]
    weights_solint.insert(main_idx, 1)
    # The required DP3 input format
    weights_solint = f'[{",".join(str(w) for w in weights_solint)}]'

    print(weights_solint)
    # return weights_solint


def main():
    get_solints(
        args.skymodel,
        solint=args.solint,
        max_flux=args.max_flux,
        min_flux=args.min_flux,
        dist=args.dist,
    )


if __name__ == "__main__":
    main()
