import logging
import re
import sys
from collections import defaultdict
from time import time
from typing import Any, Callable, Counter, Dict, List, Optional, Tuple, Union

import numpy as np
from sklearn.decomposition import PCA


def profile(func: Callable) -> Callable:
    """
    Decorator to profile a function.
    """

    def wrap(*args, **kwargs):
        fname = func.__name__
        argnames = func.__code__.co_varnames[: func.__code__.co_argcount]
        filled_args = ", ".join(
            "%s=%r" % entry
            for entry in list(zip(argnames, args[: len(argnames)]))
            + [("args", list(args[len(argnames) :]))]
            + [("kwargs", kwargs)]
        )
        logging.info(f"Started: {fname}({filled_args})")
        starting_time = time()
        output = func(*args, **kwargs)
        logging.info(f"Ended: {fname}, duration {time() - starting_time}s")
        return output

    return wrap


def scale_sample_counter(
    sample: np.ndarray, from_size: int, to_size: int
) -> Dict[str, Any]:
    """
    TODO: Add docstring
    """
    return scale_sample(sample, from_size, to_size, merge_counters)


def scale_sample(
    sample: np.ndarray,
    from_size: int,
    to_size: int,
    scaling_function: Optional[Callable] = None,
) -> Dict[str, Any]:
    """
    TODO: Add docstring
    """
    if scaling_function is None:
        logging.critical("No scaling function given!")
        sys.exit()

    if to_size is None or from_size == to_size:
        return sample

    if to_size == 0 or from_size == 0 or to_size < from_size or to_size % from_size > 0:
        logging.critical(
            f"Impossible binsize scaling requested: {from_size:,} to {to_size:,}!"
        )
        sys.exit()

    scaled_sample = dict()
    scale = to_size // from_size

    logging.info(
        f"Scaling up by a factor of {scale} from {from_size:,} to {to_size:,}."
    )

    for chrom, data in sample.items():
        new_len = int(np.ceil(len(data) / float(scale)))
        scaled_chrom = []
        for i in range(new_len):
            scaled_chrom.append(
                scaling_function(data[int(i * scale) : int((i * scale) + scale)])
            )
        scaled_sample[chrom] = np.array(scaled_chrom)

    return scaled_sample


def scale_sample_dict(
    sample: Dict[str, Any], from_size: int, to_size: int
) -> Dict[str, Any]:
    """
    TODO: Add docstring
    """
    return scale_sample(sample, from_size, to_size, merge_dicts)


def merge_counters(counters: List[Counter]) -> Counter:
    """
    TODO: add docstring
    """
    merged = Counter()
    for c in counters:
        merged.update(c)
    return merged


def merge_dicts(
    dicts: List[Dict[str, int]], defaultdict=defaultdict, int=int
) -> Dict[str, int]:
    """
    TODO add docstring
    """
    merged = defaultdict(int)
    for d in dicts:
        for k in d:
            merged[k] += d[k]
    return merged


def clip_sample(sample: Dict[str, Any], clip_lo: int = 0, clip_hi: int = 300) -> None:
    """
    TODO: Add docstring
    """

    for key in sample.keys():
        for i, counter in enumerate(sample[key]):
            sample[key][i] = clip_counter(counter, clip_lo, clip_hi)


def clip_counter(counter: Counter, lo: int, hi: int) -> Counter:
    """
    TODO: Add docstring
    """
    if lo == hi:
        return counter
    return Counter(dict(filter(lambda x: x[0] > lo and x[0] <= hi, counter.items())))


def norm_freq_mask_sample(
    sample: Dict[str, Any], cutoff: float = 0.0001, min_cutoff: int = 500
) -> None:
    """
    TODO: Add docstring
    """
    sample_freq = get_chromosome_freq(sample)
    sample_freq_full = join_chromosomes(sample_freq)
    cutoff = int(cutoff * sample_freq_full.sum())
    if min_cutoff:
        cutoff = max(cutoff, min_cutoff)
    freq_mask_sample(sample, cutoff)


def get_chromosome_freq(sample: Dict[str, Any]) -> Dict[str, np.ndarray]:
    """
    TODO: Add docstring
    """
    sample_freq = {}
    for k, v in sample.items():
        chrom_array = []
        for binv in v:
            chrom_array.append(sum(binv.values()))
        sample_freq[k] = np.array(chrom_array)
    return sample_freq


def join_chromosomes(sample: Dict[str, Any]) -> np.ndarray:
    """
    TODO: Add docstring
    """
    return np.concatenate(list(sample.values()))


def freq_mask_sample(sample: Dict[str, Any], min_observations: int) -> None:
    """
    TODO: Add docstring
    """
    n = 0
    for key in sample.keys():
        for i, counter in enumerate(sample[key]):
            n_observations = sum(counter.values())
            if n_observations >= min_observations:
                continue
            sample[key][i] = Counter()
            n += 1
    logging.info(f"Removed {n} bins that have < {min_observations} observations")


def freq_to_mean(sample: Dict[str, Any]) -> None:
    """
    TODO: Add docstring
    """
    for key in sample.keys():
        new_values = []
        for counter in sample[key]:
            mean = 0.0
            if counter:
                mean = sum(key * count for key, count in counter.items()) / sum(
                    counter.values()
                )
            new_values.append(mean)
        sample[key] = np.array(new_values)


def freq_to_median(sample: Dict[str, Any]) -> None:
    """
    TODO: Add docstring
    """
    for key in sample.keys():
        new_values = []
        for counter in sample[key]:
            median = 0.0
            if counter:
                val = np.array(list(counter.keys()))
                freq = np.array(list(counter.values()))
                ordr = np.argsort(val)
                cdf = np.cumsum(freq[ordr])
                median = val[ordr][np.searchsorted(cdf, cdf[-1] // 2)]
            new_values.append(median)
        sample[key] = np.array(new_values)


def convert_txt(text: str) -> Union[int, str]:
    """
    Returns an integer if the string is a number, otherwise returns the string in lowercase.
    """
    return int(text) if text.isdigit() else text.lower()


def natural_sort(keys_list: List[str]) -> List[str]:
    """
    Sorts the given iterable in a natural way.
    """
    return sorted(
        keys_list, key=lambda key: [convert_txt(c) for c in re.split("([0-9]+)", key)]
    )


def get_mask(samples: np.ndarray, rc: bool = False) -> Tuple[np.ndarray, List[int]]:
    """
    TODO: Add docstring
    """
    by_chr = []
    bins_per_chr = []
    sample_count = len(samples)

    for chrom in natural_sort(samples[0].keys()):
        max_len = max([sample[chrom].shape[0] for sample in samples])
        this_chr = np.zeros((max_len, sample_count), dtype=float)
        bins_per_chr.append(max_len)

        for i, sample in enumerate(samples):
            this_chr[:, i] = sample[chrom]
        by_chr.append(this_chr)

    all_data = np.concatenate(by_chr, axis=0)

    if rc:
        sum_per_sample = np.sum(all_data, 0)
        all_data = all_data / sum_per_sample

    sum_per_bin = np.sum(all_data, 1)
    mask = sum_per_bin > 0

    return mask, bins_per_chr


def train_pca(ref_data: np.ndarray, pcacomp: int = 1) -> Tuple[np.ndarray, PCA]:
    """
    TODO: Add docstring
    """
    t_data = ref_data.T
    pca = PCA(n_components=pcacomp)
    pca.fit(t_data)
    PCA(copy=True, whiten=False)
    transformed = pca.transform(t_data)
    inversed = pca.inverse_transform(transformed)
    corrected = t_data / inversed
    return corrected.T, pca


def normalize_and_mask(
    samples: List[Dict[str, Any]], chrs: Union[range, List[int]], mask, rc: bool = False
) -> np.ndarray:
    """
    TODO: Add docstring
    """
    by_chr = []
    sample_count = len(samples)

    for chrom in chrs:
        max_len = max([sample[str(chrom)].shape[0] for sample in samples])
        this_chr = np.zeros((max_len, sample_count), dtype=float)
        for i, sample in enumerate(samples):
            this_chr[:, i] = sample[str(chrom)]
        by_chr.append(this_chr)

    all_data = np.concatenate(by_chr, axis=0)

    if rc:
        sum_per_sample = np.sum(all_data, 0)
        all_data = all_data / sum_per_sample

    masked_data = all_data[mask, :]

    return masked_data


def split_by_chr(start, end, chr_bin_sums):
    areas = []
    tmp = [0, start, 0]
    for i, val in enumerate(chr_bin_sums):
        tmp[0] = i
        if val >= end:
            break
        if start < val < end:
            tmp[2] = val
            areas.append(tmp)
            tmp = [i, val, 0]
        tmp[1] = val
    tmp[2] = end
    areas.append(tmp)
    return areas
