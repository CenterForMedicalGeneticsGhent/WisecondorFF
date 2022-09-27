import logging
import os
import warnings
from pathlib import Path
from typing import List, Optional

import numpy as np
import typer

from wisecondorff.convert import convert_reads
from wisecondorff.detect import (
    detect_wrap,
    exec_cbs,
    generate_segments,
    get_res_to_nparray,
    res_to_nestedlist,
)
from wisecondorff.reference import reference_construct, reference_prep
from wisecondorff.utils import (
    clip_sample,
    freq_to_mean,
    get_mask,
    norm_freq_mask_sample,
    profile,
    scale_sample,
    scale_sample_counter,
)

# from wisecondorff import __version__

__author__ = "Tom Okveld, Matthias De Smet"
__copyright__ = "Tom Okveld, Center for Medical Genetics, Ghent University, Belgium"
__license__ = "Attribution-NonCommercial-ShareAlike 4.0 International"

warnings.filterwarnings("ignore")

app = typer.Typer()


@profile
@app.command()
def convert(
    infile: Path = typer.Argument(
        ..., help="Input BAM/CRAM file for conversion.", show_default=False
    ),
    outfile: Path = typer.Argument(..., help="Output .npz file.", show_default=False),
    binsize: int = typer.Option(default=50000, help="Bin size for conversion."),
    min_mapq: int = typer.Option(default=1, help="Minimum mapping quality."),
    reference_fasta: Optional[Path] = typer.Option(
        default=None, help="Reference fasta for CRAM conversion."
    ),
) -> None:

    """
    Read and process a BAM or CRAM file.
    """
    sample, qual_info = convert_reads(infile, binsize, min_mapq, reference_fasta)
    np.savez_compressed(
        outfile,
        quality=qual_info,
        args={
            "binsize": binsize,
            "map_quality": min_mapq,
            "infile": str(infile),
            "outfile": str(outfile),
        },
        sample=sample,
    )


@profile
@app.command()
def reference(
    npz_list: List[Path] = typer.Argument(
        ..., help="List of npz files", show_default=False
    ),
    out: Path = typer.Argument(
        ...,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
        show_default=False,
    ),
    refsize: int = typer.Option(
        default=300, help="Number of reference regions per region"
    ),
    binsize: int = typer.Option(
        default=500000,
        help="Scale samples to this region size (multiples of existing region size only)",
    ),
    fs_clip_low: int = typer.Option(
        default=0,
        help="Lower bound for the inclusion range of the fragment size distribution",
    ),
    fs_clip_high: int = typer.Option(
        default=3000,
        help="Upper bound for the inclusion range of the fragment size distribution",
    ),
    rc_clip_norm: float = typer.Option(
        default=0.0001,
        help="Lower bound cutoff based on the normalized read count across regions (regions that fall below this cutoff are masked)",
    ),
    rc_clip_abs: int = typer.Option(
        default=500,
        help="Lower bound cutoff based on the absolute read count across regions (regions that fall below this cutoff are masked)",
    ),
) -> None:
    """
    Construct a reference from a list of npz files
    """
    split_path = list(os.path.split(out))
    if split_path[-1][-4:] == ".npz":
        split_path[-1] = split_path[-1][:-4]

    rc_samples = []
    fs_samples = []

    for npz in npz_list:
        logging.info(f"Loading: {npz}")
        npzdata = np.load(npz, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        sample_args = npzdata["args"].item()
        sample_binsize = int(sample_args["binsize"])

        fs_sample = scale_sample_counter(sample["FS"], sample_binsize, binsize)
        clip_sample(fs_sample, clip_lo=fs_clip_low, clip_hi=fs_clip_high)
        norm_freq_mask_sample(fs_sample, cutoff=rc_clip_norm, min_cutoff=rc_clip_abs)
        freq_to_mean(fs_sample)

        fs_samples.append(fs_sample)
        rc_samples.append(scale_sample(sample["RC"], sample_binsize, binsize, np.sum))

    fs_samples = np.array(fs_samples)
    rc_samples = np.array(rc_samples)

    fs_total_mask, fs_bins_per_chr = get_mask(fs_samples, rc=False)
    rc_total_mask, rc_bins_per_chr = get_mask(rc_samples, rc=True)

    fs_auto = reference_prep(
        binsize, refsize, fs_samples, fs_total_mask, fs_bins_per_chr, rc=False
    )
    rc_auto = reference_prep(
        binsize, refsize, rc_samples, rc_total_mask, rc_bins_per_chr, rc=True
    )

    reference_construct(fs_auto, refsize)
    reference_construct(rc_auto, refsize)

    np.savez_compressed(out, reference={"RC": rc_auto, "FS": fs_auto})


@profile
@app.command()
def detect(
    in_npz: Path = typer.Argument(
        ..., help="Input sample npz file", show_default=False
    ),
    reference_npz: Path = typer.Argument(
        ..., help="Reference npz file", show_default=False
    ),
    minrefbins: int = typer.Option(
        default=150, help="Minimum amount of sensible reference bins per target bin."
    ),
    maskrepeats: int = typer.Option(
        default=5,
        help="Regions with distances > mean + sd * 3 will be masked. Number of masking cycles.",
    ),
    zscore: float = typer.Option(
        default=5, help="Z-score cut-off for aberration calling."
    ),
    fs_clip_low: int = typer.Option(
        default=0,
        help="Lower bound for the inclusion range of the fragment size distribution",
    ),
    fs_clip_high: int = typer.Option(
        default=3000,
        help="Upper bound for the inclusion range of the fragment size distribution",
    ),
    rc_clip_norm: float = typer.Option(
        default=0.0001,
        help="Lower bound cutoff based on the normalized read count across regions (regions that fall below this cutoff are masked)",
    ),
    rc_clip_abs: int = typer.Option(
        default=500,
        help="Lower bound cutoff based on the absolute read count across regions (regions that fall below this cutoff are masked)",
    ),
) -> None:
    """
    Detect copy number aberrations in a sample using a reference.
    """
    ref_npz = np.load(str(reference_npz), encoding="latin1", allow_pickle=True)
    sample_npz = np.load(str(in_npz), encoding="latin1", allow_pickle=True)

    ref = ref_npz["reference"].item()

    sample_args = sample_npz["args"].item()
    sample_data = sample_npz["sample"].item()

    rc_sample = sample_data["RC"]
    fs_sample = sample_data["FS"]

    rc_sample = scale_sample(
        rc_sample, sample_args["binsize"], ref["RC"]["binsize"], np.sum
    )
    fs_sample = scale_sample_counter(
        fs_sample, sample_args["binsize"], ref["FS"]["binsize"]
    )

    clip_sample(fs_sample, clip_lo=fs_clip_low, clip_hi=fs_clip_high)
    norm_freq_mask_sample(fs_sample, cutoff=rc_clip_norm, min_cutoff=rc_clip_abs)
    freq_to_mean(fs_sample)

    rc_results, _ = detect_wrap(
        rc_sample,
        ref["RC"],
        rc=True,
        maskrepeats=maskrepeats,
        minrefbins=minrefbins,
        zscore=zscore,
    )
    fs_results, fs_rem_input = detect_wrap(
        fs_sample,
        ref["FS"],
        rc=False,
        maskrepeats=maskrepeats,
        minrefbins=minrefbins,
        zscore=zscore,
    )

    _rc_results = get_res_to_nparray(rc_results)
    _fs_results = get_res_to_nparray(fs_results)

    comb_results = {}
    comb_results["results_z"] = (
        _rc_results["results_z"] + (-_fs_results["results_z"])
    ) / np.sqrt(2)
    comb_results["results_r"] = (
        _rc_results["results_r"] + (-_fs_results["results_r"])
    ) / np.sqrt(2)
    comb_results["results_w"] = (
        _rc_results["results_w"] + (_fs_results["results_w"])
    ) / 2
    comb_results["results_nr"] = (
        _rc_results["results_nr"] + (_fs_results["results_nr"])
    ) / 2

    res_to_nestedlist(comb_results)

    comb_results["results_c"] = exec_cbs(fs_rem_input, comb_results)

    generate_segments(fs_rem_input, comb_results)


######
# def make_generic(cvalue, func_xs, comp_xs, typecast, message):
#     """Function generator that can handle most of the different types needed by argparse

#     Args:
#         cvalue (T): Optional value to compare to with the argument of the inner function
#         func_xs (List[func(x: T) -> bool]): Optional list value of functions that are applied onto the argument of the inner function
#         comp_xs (List[func(x: T, y: T) -> bool] or List[func(x: T) -> bool]):
#         List value of functions (two-argument if cvalue is set, otherwise only applied on the argument of the inner function)
#         typecast (T): Typecast argument to cast the argument of the inner function
#         message (str): fstring used in error messages should follow the format of {tvalue} {cvalue}
#     Returns:
#         f(x: T) -> x: T
#     """

#     def check_generic(value):
#         if func_xs:
#             tvalue = functools.reduce(lambda res, f: f(res), func_xs, typecast(value))
#         else:
#             tvalue = typecast(value)

#         if cvalue is None:
#             if functools.reduce(lambda res, f: f(res), comp_xs, tvalue):
#                 raise argparse.ArgumentTypeError(message.format(tvalue, cvalue))
#         else:
#             if any([func(tvalue, cvalue) for func in comp_xs]):
#                 raise argparse.ArgumentTypeError(message.format(tvalue, cvalue))

#         return typecast(value)

#     return check_generic


# REGION_MIN = 5000
# REGION_MAX = 20000000

# QUAL_MIN = 0
# QUAL_MAX = 30

# check_bins = make_generic(
#     REGION_MAX,
#     None,
#     [operator.ge, lambda *x: x[0] < REGION_MIN],
#     int,
#     '"{}" invalid region size set, should be between 5000 bp to {} MB',
# )

# check_qual = make_generic(
#     QUAL_MAX,
#     None,
#     [operator.ge, lambda *x: x[0] < QUAL_MIN],
#     int,
#     '"{}" invalid mapping quality set, should be between 0 bp to 30',
# )


if __name__ == "__main__":
    app()
