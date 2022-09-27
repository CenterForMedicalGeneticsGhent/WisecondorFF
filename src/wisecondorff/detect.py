import json
import logging
import os
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.stats import norm
from sklearn.decomposition import PCA


def generate_segments(rem_input, results) -> None:
    for segment in results["results_c"]:
        chr_name = str(segment[0] + 1)
        row = [
            chr_name,
            int(segment[1] * rem_input["binsize"] + 1),
            int(segment[2] * rem_input["binsize"]),
            segment[4],
            segment[3],
        ]
        if float(segment[3]) > rem_input["zscore"]:
            print("{}\tgain\n".format("\t".join([str(x) for x in row])))
        elif float(segment[3]) < -rem_input["zscore"]:
            print("{}\tloss\n".format("\t".join([str(x) for x in row])))


def res_to_nestedlist(results) -> None:
    for k, v in results.items():
        results[k] = [list(i) for i in v]


def exec_R(json_dict) -> Dict[str, Any]:
    json.dump(json_dict, open(json_dict["infile"], "w"))

    r_cmd = ["Rscript", json_dict["R_script"], "--infile", json_dict["infile"]]
    logging.debug(f"CBS cmd: {r_cmd}")

    try:
        subprocess.check_call(r_cmd)
    except subprocess.CalledProcessError as e:
        logging.critical(f"Rscript failed: {e}")
        sys.exit()

    os.remove(json_dict["infile"])
    if "outfile" in json_dict.keys():
        json_out = json.load(open(json_dict["outfile"]))
        os.remove(json_dict["outfile"])
        return json_out


def exec_cbs(rem_input, results):
    json_cbs_dir = os.path.abspath("test_" + str(uuid.uuid4()) + "_CBS_tmp")

    json_dict = {
        "R_script": str(
            Path(os.path.dirname(os.path.realpath(__file__))).joinpath("include/CBS.R")
        ),
        "ref_gender": "A",
        "alpha": str(1e-4),
        "binsize": str(rem_input["binsize"]),
        "results_r": results["results_r"],
        "results_w": results["results_w"],
        "infile": str(f"{json_cbs_dir}_01.json"),
        "outfile": str(f"{json_cbs_dir}_02.json"),
    }

    results_c = get_processed_cbs(exec_R(json_dict))
    segment_z = get_z_score(results_c, results)
    results_c = [
        results_c[i][:3] + [segment_z[i]] + [results_c[i][3]]
        for i in range(len(results_c))
    ]
    return results_c


def get_processed_cbs(cbs_data):
    results_c = []
    for i, segment in enumerate(cbs_data):
        chr = int(segment["chr"]) - 1
        s = int(segment["s"])
        e = int(segment["e"])
        r = segment["r"]
        results_c.append([chr, s, e, r])
    return results_c


def get_z_score(results_c, results):
    results_nr, results_r, results_w = (
        results["results_nr"],
        results["results_r"],
        results["results_w"],
    )
    zs = []
    for segment in results_c:
        segment_nr = results_nr[segment[0]][segment[1] : segment[2]]
        segment_rr = results_r[segment[0]][segment[1] : segment[2]]
        segment_nr = [
            segment_nr[i] for i in range(len(segment_nr)) if segment_rr[i] != 0
        ]
        segment_w = results_w[segment[0]][segment[1] : segment[2]]
        segment_w = [segment_w[i] for i in range(len(segment_w)) if segment_rr[i] != 0]
        null_segments = [
            np.ma.average(x, weights=segment_w) for x in np.transpose(segment_nr)
        ]
        null_mean = np.ma.mean([x for x in null_segments if np.isfinite(x)])
        null_sd = np.ma.std([x for x in null_segments if np.isfinite(x)])
        z = (segment[3] - null_mean) / null_sd
        z = min(z, 1000)
        z = max(z, -1000)
        zs.append(z)
    return zs


def detect_wrap(
    sample,
    final_dict,
    rc: bool = False,
    maskrepeats: int = 5,
    minrefbins: int = 150,
    zscore: float = 5,
):
    results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
        maskrepeats, sample, final_dict, rc=rc
    )

    null_ratios_aut_per_bin = final_dict["null_ratios"]

    rem_input = {
        "binsize": int(final_dict["binsize"]),
        "zscore": zscore,
        "mask": final_dict["mask"],
        "bins_per_chr": np.array(final_dict["bins_per_chr"]),
        "masked_bins_per_chr": np.array(final_dict["masked_bins_per_chr"]),
        "masked_bins_per_chr_cum": np.array(final_dict["masked_bins_per_chr_cum"]),
    }

    results_z = results_z - m_z
    results_w = results_w / np.nanmean(results_w)

    if np.isnan(results_w).any() or np.isinf(results_w).any():
        logging.warning(
            "Non-numeric values found in weights -- reference too small."
            "Circular binary segmentation and z-scoring will be unweighted"
        )
        results_w = np.ones(len(results_w))

    null_ratios = np.array([x.tolist() for x in null_ratios_aut_per_bin])

    results = {
        "results_r": results_r,
        "results_z": results_z,
        "results_w": results_w,
        "results_nr": null_ratios,
        "ref_sizes": ref_sizes,
    }

    for result in results.keys():
        results[result] = get_post_processed_result(
            minrefbins, results[result], ref_sizes, rem_input
        )

    log_trans(results, m_lr)

    return results, rem_input


def normalize(maskrepeats, sample, final_dict, rc=False):
    sample = coverage_normalize_and_mask(sample, final_dict, rc)
    sample = project_pc(sample, final_dict)
    results_w = get_weights(final_dict)
    optimal_cutoff = get_optimal_cutoff(final_dict, maskrepeats)
    results_z, results_r, ref_sizes, m_lr, m_z = normalize_repeat(
        sample, final_dict, optimal_cutoff
    )
    return results_r, results_z, results_w, ref_sizes, m_lr, m_z


def get_weights(ref_file):
    inverse_weights = [np.mean(np.sqrt(x)) for x in ref_file["distances"]]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


def get_optimal_cutoff(ref_file, repeats):
    distances = ref_file["distances"]
    cutoff = float("inf")
    for i in range(0, repeats):
        mask = distances < cutoff
        average = np.average(distances[mask])
        stddev = np.std(distances[mask])
        cutoff = average + 3 * stddev
    return cutoff


def coverage_normalize_and_mask(sample, ref_file, rc=False):
    by_chr = []
    chromosomes = range(1, len(ref_file["bins_per_chr"]) + 1)

    for chrom in chromosomes:
        this_chr = np.zeros(ref_file["bins_per_chr"][chrom - 1], dtype=float)
        min_len = min(ref_file["bins_per_chr"][chrom - 1], len(sample[str(chrom)]))
        this_chr[:min_len] = sample[str(chrom)][:min_len]
        by_chr.append(this_chr)
    all_data = np.concatenate(by_chr, axis=0)

    if rc:
        all_data = all_data / np.sum(all_data)

    masked_data = all_data[ref_file["mask"]]

    return masked_data


def project_pc(sample_data, ref_file):
    pca = PCA(n_components=ref_file["pca_components"].shape[0])
    pca.components_ = ref_file["pca_components"]
    pca.mean_ = ref_file["pca_mean"]

    transform = pca.transform(np.array([sample_data]))

    reconstructed = np.dot(transform, pca.components_) + pca.mean_
    reconstructed = reconstructed[0]
    return sample_data / reconstructed


def normalize_repeat(
    test_data, ref_file, optimal_cutoff
) -> Tuple[Any, Any, Any, Any, Any]:
    test_copy = np.copy(test_data)
    for i in range(3):  # TODO: Parameter
        results_z, results_r, ref_sizes = _normalize_once(
            test_data, test_copy, ref_file, optimal_cutoff
        )

        test_copy[np.abs(results_z) >= norm.ppf(0.99)] = -1
    m_lr = np.nanmedian(np.log2(results_r))
    m_z = np.nanmedian(results_z)

    return results_z, results_r, ref_sizes, m_lr, m_z


def _normalize_once(test_data, test_copy, ref_file, optimal_cutoff):
    masked_bins_per_chr = ref_file["masked_bins_per_chr"]
    masked_bins_per_chr_cum = ref_file["masked_bins_per_chr_cum"]
    results_z = np.zeros(masked_bins_per_chr_cum[-1])
    results_r = np.zeros(masked_bins_per_chr_cum[-1])
    ref_sizes = np.zeros(masked_bins_per_chr_cum[-1])
    indexes = ref_file["indexes"]
    distances = ref_file["distances"]

    i = 0
    i2 = 0
    for chrom in list(range(len(masked_bins_per_chr))):
        start = masked_bins_per_chr_cum[chrom] - masked_bins_per_chr[chrom]
        end = masked_bins_per_chr_cum[chrom]
        chr_data = np.concatenate(
            (
                test_copy[
                    : masked_bins_per_chr_cum[chrom] - masked_bins_per_chr[chrom]
                ],
                test_copy[masked_bins_per_chr_cum[chrom] :],
            )
        )
        for index in indexes[start:end]:
            ref_data = chr_data[index[distances[i] < optimal_cutoff]]
            ref_data = ref_data[ref_data >= 0]
            ref_stdev = np.std(ref_data)

            results_z[i2] = (test_data[i] - np.mean(ref_data)) / ref_stdev
            results_r[i2] = test_data[i] / np.median(ref_data)
            ref_sizes[i2] = ref_data.shape[0]
            i += 1
            i2 += 1

    return results_z, results_r, ref_sizes


def get_post_processed_result(minrefbins, result, ref_sizes, rem_input) -> List[Any]:
    infinite_mask = ref_sizes < minrefbins
    result[infinite_mask] = 0
    inflated_results = inflate_results(result, rem_input)

    final_results = []
    for chr in range(len(rem_input["bins_per_chr"])):
        chr_data = inflated_results[
            sum(rem_input["bins_per_chr"][:chr]) : sum(
                rem_input["bins_per_chr"][: chr + 1]
            )
        ]
        final_results.append(chr_data)

    return final_results


def inflate_results(results, rem_input):
    temp = [0 for x in rem_input["mask"]]
    j = 0
    for i, val in enumerate(rem_input["mask"]):
        if val:
            temp[i] = results[j]
            j += 1
    return temp


def get_res_to_nparray(results):
    new_results = {}
    new_results["results_z"] = np.array([np.array(x) for x in results["results_z"]])
    new_results["results_w"] = np.array([np.array(x) for x in results["results_w"]])
    new_results["results_nr"] = np.array([np.array(x) for x in results["results_nr"]])
    new_results["results_r"] = np.array([np.array(x) for x in results["results_r"]])
    return new_results


def log_trans(results, log_r_median) -> None:
    for chr in range(len(results["results_r"])):
        results["results_r"][chr] = np.log2(results["results_r"][chr])

    results["results_r"] = [x.tolist() for x in results["results_r"]]

    for c in range(len(results["results_r"])):
        for i, rR in enumerate(results["results_r"][c]):
            if not np.isfinite(rR):
                results["results_r"][c][i] = 0
                results["results_z"][c][i] = 0
                results["results_w"][c][i] = 0
            if results["results_r"][c][i] != 0:
                results["results_r"][c][i] = results["results_r"][c][i] - log_r_median
