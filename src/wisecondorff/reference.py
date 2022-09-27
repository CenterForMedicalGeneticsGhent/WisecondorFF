import random
from typing import Any, Dict, Tuple

import numpy as np

from wisecondorff.utils import normalize_and_mask, split_by_chr, train_pca


def reference_prep(
    binsize: int, refsize: int, samples, mask, bins_per_chr, rc: bool = False
) -> Dict[str, Any]:
    """
    TODO: Add docstring
    """
    bins_per_chr = bins_per_chr[:22]
    mask = mask[: np.sum(bins_per_chr)]

    masked_data = normalize_and_mask(samples, range(1, 23), mask, rc)
    pca_corrected_data, pca = train_pca(masked_data)

    masked_bins_per_chr = [
        sum(mask[sum(bins_per_chr[:i]) : sum(bins_per_chr[:i]) + x])
        for i, x in enumerate(bins_per_chr)
    ]
    masked_bins_per_chr_cum = [
        sum(masked_bins_per_chr[: x + 1]) for x in range(len(masked_bins_per_chr))
    ]

    return {
        "binsize": binsize,
        "refsize": refsize,
        "mask": mask,
        "masked_data": masked_data,
        "bins_per_chr": bins_per_chr,
        "masked_bins_per_chr": masked_bins_per_chr,
        "masked_bins_per_chr_cum": masked_bins_per_chr_cum,
        "pca_corrected_data": pca_corrected_data,
        "pca_components": pca.components_,  # type: ignore
        "pca_mean": pca.mean_,  # type: ignore
    }


def get_reference(
    pca_corrected_data: np.ndarray,
    masked_bins_per_chr: Dict[str, int],
    masked_bins_per_chr_cum,
    ref_size: int,
) -> Tuple[Any, Any, Any]:
    big_indexes = []
    big_distances = []

    regions = split_by_chr(0, masked_bins_per_chr_cum[-1], masked_bins_per_chr_cum)
    for (chrom, start, end) in regions:
        chr_data = np.concatenate(
            (
                pca_corrected_data[
                    : masked_bins_per_chr_cum[chrom] - masked_bins_per_chr[chrom], :
                ],
                pca_corrected_data[masked_bins_per_chr_cum[chrom] :, :],
            )
        )

        part_indexes, part_distances = get_ref_for_bins(
            ref_size, start, end, pca_corrected_data, chr_data
        )

        big_indexes.extend(part_indexes)
        big_distances.extend(part_distances)

    index_array = np.array(big_indexes)
    distance_array = np.array(big_distances)
    null_ratio_array = np.zeros(
        (len(distance_array), min(len(pca_corrected_data[0]), 100))
    )  # TODO: make parameter
    samples = np.transpose(pca_corrected_data)

    for null_i, case_i in enumerate(
        random.sample(
            range(len(pca_corrected_data[0])), min(len(pca_corrected_data[0]), 100)
        )
    ):
        sample = samples[case_i]
        for bin_i in list(range(len(sample))):
            r = np.log2(sample[bin_i] / np.median(sample[index_array[bin_i]]))
            null_ratio_array[bin_i][null_i] = r

    return index_array, distance_array, null_ratio_array


def get_ref_for_bins(ref_size, start, end, pca_corrected_data, chr_data):
    ref_indexes = np.zeros((end - start, ref_size), dtype=np.int32)
    ref_distances = np.ones((end - start, ref_size))
    for cur_bin in range(start, end):
        bin_distances = np.sum(
            np.power(chr_data - pca_corrected_data[cur_bin, :], 2), 1
        )

        unsrt_ranked_idx = np.argpartition(bin_distances, ref_size)[:ref_size]
        ranked_idx = unsrt_ranked_idx[np.argsort(bin_distances[unsrt_ranked_idx])]
        ranked_distances = bin_distances[ranked_idx]

        ref_indexes[cur_bin - start, :] = ranked_idx
        ref_distances[cur_bin - start, :] = ranked_distances

    return ref_indexes, ref_distances


def reference_construct(ref_dict: Dict[str, Any], ref_size: int = 300) -> None:
    """
    TODO: Add docstring
    """
    indexes, distances, null_ratios = get_reference(
        ref_dict["pca_corrected_data"],
        ref_dict["masked_bins_per_chr"],
        ref_dict["masked_bins_per_chr_cum"],
        ref_size=ref_size,
    )
    ref_dict["indexes"] = indexes
    ref_dict["distances"] = distances
    ref_dict["null_ratios"] = null_ratios
