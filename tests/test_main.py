# import pytest
from typer.testing import CliRunner

from wisecondorff.main import (
    app,
    clip_counter,
    clip_sample,
    convert,
    convert_reads,
    convert_txt,
    coverage_normalize_and_mask,
    detect,
    exec_cbs,
    exec_R,
    freq_mask_sample,
    freq_to_mean,
    freq_to_median,
    generate_segments,
    get_chromosome_freq,
    get_mask,
    get_optimal_cutoff,
    get_post_processed_result,
    get_ref_for_bins,
    get_reference,
    get_res_to_nparray,
    get_weights,
    get_z_score,
    inflate_results,
    join_chromosomes,
    log_trans,
    make_generic,
    merge_counters,
    merge_dicts,
    natural_sort,
    norm_freq_mask_sample,
    normalize,
    normalize_and_mask,
    normalize_repeat,
    profile,
    project_pc,
    read_pair_gen,
    reference,
    reference_construct,
    reference_prep,
    res_to_nestedlist,
    scale_sample,
    scale_sample_array,
    scale_sample_counter,
    scale_sample_dict,
    signal_kill,
    split_by_chr,
    train_pca,
)

__author__ = "Tom Okveld, Matthias De Smet"
__copyright__ = "Tom Okveld, Center for Medical Genetics, Ghent University, Belgium"
__license__ = "Attribution-NonCommercial-ShareAlike 4.0 International"


runner = CliRunner()


def test_app():
    """
    Test the main CLI.
    """
    res = runner.invoke(app, ["--help"])
    assert res.exit_code == 0


def test_app_convert():
    """
    Test the convert CLI.
    """
    res = runner.invoke(app, ["convert", "--help"])
    assert res.exit_code == 0


def test_app_reference():
    """
    Test the reference CLI.
    """
    res = runner.invoke(app, ["reference", "--help"])
    assert res.exit_code == 0


def test_app_detect():
    """
    Test the detect CLI.
    """
    res = runner.invoke(app, ["detect", "--help"])
    assert res.exit_code == 0


def test_convert():
    """
    Test the main convert function.
    """
    convert()
    pass


def test_reference():
    """
    Test the main reference function
    """
    reference()
    pass


def test_detect():
    """
    Test the main detect function
    """
    detect()
    pass


def test_profile():
    profile()
    pass  # TODO: test profile


def test_signal_kill():
    signal_kill()
    pass  # TODO: test signal_kill


def test_read_pair_gen():
    read_pair_gen()
    pass  # TODO: test read_pair_gen


def test_convert_reads():
    convert_reads()
    pass  # TODO: test convert_reads


def test_convert_txt():
    convert_txt()
    pass  # TODO: test wcr_convert


def test_scale_sample():
    scale_sample()
    pass  # TODO: test scale_sample


def test_scale_sample_array():
    scale_sample_array()
    pass  # TODO: test scale_sample_array


def test_merge_counters():
    merge_counters()
    pass  # TODO: test merge_counters


def test_merge_dicts():
    merge_dicts()
    pass  # TODO: test merge_dicts


def test_scale_sample_counter():
    scale_sample_counter()
    pass  # TODO: test scale_sample_counter


def test_scale_sample_dict():
    scale_sample_dict()
    pass  # TODO: test scale_sample_dict


def test_clip_counter():
    clip_counter()
    pass  # TODO: test clip_counter


def test_clip_sample():
    clip_sample()
    pass  # TODO: test clip_sample


def test_get_chromosome_freq():
    get_chromosome_freq()
    pass  # TODO: test get_chromosome_freq


def test_join_chromosomes():
    join_chromosomes()
    pass  # TODO: test join_chromosomes


def test_freq_mask_sample():
    freq_mask_sample()
    pass  # TODO: test freq_mask_sample


def test_norm_freq_mask_sample():
    norm_freq_mask_sample()
    pass  # TODO: test norm_freq_mask_sample


def test_freq_to_mean():
    freq_to_mean()
    pass  # TODO: test freq_to_mean


def test_freq_to_median():
    freq_to_median()
    pass  # TODO: test freq_to_median


def test_natural_sort():
    natural_sort()
    pass  # TODO: test natural_sort


def test_get_mask():
    get_mask()
    pass  # TODO: test get_mask


def test_train_pca():
    train_pca()
    pass  # TODO: test train_pca


def test_normalize_and_mask():
    normalize_and_mask()
    pass  # TODO: test normalize_and_mask


def test_reference_prep():
    reference_prep()
    pass  # TODO: test reference_prep


def test_get_reference():
    get_reference()
    pass  # TODO: test get_reference


def test_get_ref_for_bins():
    get_ref_for_bins()
    pass  # TODO: test get_ref_for_bins


def test_split_by_chr():
    split_by_chr()
    pass  # TODO: test split_by_chr


def test_reference_construct():
    reference_construct()
    pass  # TODO: test reference_construct


def test_coverage_normalize_and_mask():
    coverage_normalize_and_mask()
    pass  # TODO: test coverage_normalize_and_mask


def test_project_pc():
    project_pc()
    pass  # TODO: test project_pc


def test_get_weights():
    get_weights()
    pass  # TODO: test get_weights


def test_get_optimal_cutoff():
    get_optimal_cutoff()
    pass  # TODO: test get_optimal_cutoff


def test_normalize_repeat():
    normalize_repeat()
    pass  # TODO: test normalize_repeat


def test_normalize():
    normalize()
    pass  # TODO: test normalize


def test_inflate_results():
    inflate_results()
    pass  # TODO: test inflate_results


def test_get_post_processed_result():
    get_post_processed_result()
    pass  # TODO: test get_post_processed_result


def test_log_trans():
    log_trans()
    pass  # TODO: test log_trans


def test_get_res_to_nparray():
    get_res_to_nparray()
    pass  # TODO: test get_res_to_nparray


def test_res_to_nestedlist():
    res_to_nestedlist()
    pass  # TODO: test res_to_nestedlist


def test_get_z_score():
    get_z_score()
    pass  # TODO: test get_z_score


def test_exec_R():
    exec_R()
    pass  # TODO: test exec_R


def test_exec_cbs():
    exec_cbs()
    pass  # TODO: test exec_cbs


def test_generate_segments():
    generate_segments()
    pass  # TODO: test generate_segments


def test_make_generic():
    make_generic()
    pass  # TODO: test make_generic
