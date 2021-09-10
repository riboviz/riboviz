"""
Integration test suite.

The integration test suite runs
:py:const:`riboviz.test.NEXTFLOW_WORKFLOW` (via Nextflow) using a
given configuration file, then compares the results to a directory of
pre-calculated results, specified by the user.

Usage::

    pytest riboviz/test/integration/test_integration.py
      --expected=DIRECTORY
      [--skip-workflow]
      [--check-index-tmp]
      [--config-file=FILE]

See :py:mod:`riboviz.test.integration.conftest` for information on the
command-line parameters and the fixtures used by these tests.

If the configuration uses environment variable tokens, then these
should be defined in the bash shell within which ``pytest`` is
run. Alternatively, they can be provided when running ``pytest``, for
example::

    $ RIBOVIZ_SAMPLES=data/ RIBOVIZ_ORGANISMS=vignette/input/ \
      RIBOVIZ_DATA=data/ pytest ...

As the expected data directories and those with the data to be tested
may vary in their paths the following approach is used:

* The paths of the directories with the data to be tested are taken to
  be those specified in the configuration file.
* The paths of the directories with the expected data are taken to be
  relative to the ``--expected`` directory and to share common names
  with the final directories of each path of the actual data
  directories.

For example, if the configuration file has::

    dir_index: vignette/index
    dir_out: vignette/simdata_umi_output
    dir_tmp: vignette/simdata_umi_tmp

and ``--expected`` is ``/home/user/simdata-umi-data`` then directories
with the data to be tested are::

    vignette/index
    vignette/simdata_umi_output
    vignette/simdata_umi_tmp

and the directories with the expected data are::

    /home/user/simdata-umi-data/index
    /home/user/simdata-umi-data/simdata_umi_output
    /home/user/simdata-umi-data/simdata_umi_tmp

If running with a configuration that used UMI extraction,
deduplication and grouping then note that:

* UMI deduplication statistics files (files prefixed
  by :py:const:`riboviz.workflow_files.DEDUP_STATS_PREFIX`) can differ
  between runs depending on which reads are removed by ``umi_tools
  dedup``, so only the existence of the files is checked.
* UMI group file post-deduplication files,
  (:py:const:`riboviz.workflow_files.POST_DEDUP_GROUPS_TSV`) can
  differ between runs depending on which reads are removed by
  ``umi_tools dedup``, so only the existence of the file is checked.
* BAM file output by deduplication (``<SAMPLE>.bam``) can differ
  between runs depending on which reads are removed by ``umi_tools
  dedup``, so only the existence of the file is checked.
"""
import os
import pytest
import pysam
from riboviz import bedgraph
from riboviz import count_reads as count_reads_module
from riboviz import environment
from riboviz import fastq
from riboviz import h5
from riboviz import html
from riboviz import hisat2
from riboviz import sam_bam
from riboviz import utils
from riboviz import workflow_files
from riboviz import workflow_r
from riboviz.test import nextflow
from riboviz import test


@pytest.fixture(scope="module")
def prep_riboviz_fixture(skip_workflow_fixture, config_fixture):
    """
    Run :py:const:`riboviz.test.NEXTFLOW_WORKFLOW` (via Nextflow)
    if ``skip_workflow_fixture`` is not ``True``.

    :param skip_workflow_fixture: Should workflow not be run?
    :type skip_workflow_fixture: bool
    :param config_fixture: Configuration file
    :type config_fixture: str or unicode
    """
    if not skip_workflow_fixture:
        env_vars = environment.get_environment_vars()
        exit_code = nextflow.run_nextflow(config_fixture,
                                          envs=env_vars)
        assert exit_code == 0, \
            "prep_riboviz returned non-zero exit code %d" % exit_code


@pytest.fixture(scope="function")
def scratch_directory(tmpdir):
    """
    Create a scratch directory.

    :param tmpdir: Temporary directory (pytest built-in fixture)
    :type tmpdir: py._path.local.LocalPath
    :return: directory
    :rtype: py._path.local.LocalPath
    """
    return tmpdir.mkdir("scratch")


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("index", list(range(1, test.NUM_INDICES)))
def test_hisat2_build_index(build_indices, expected_fixture, dir_index,
                            index_prefix, index):
    """
    Test ``hisat2-build`` index file sizes for equality. See
    :py:func:`riboviz.utils.equal_file_sizes`.

    Skipped if :py:const:`riboviz.params.BUILD_INDICES` is ``False``.

    :param build_indicex: Configuration parameter
    :type build_indices: boolean
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_index: Index files directory
    :type dir_index: str or unicode
    :param index_prefix: Index file name prefix
    :type index_prefix: str or unicode
    :param index: File name index
    :type index: int
    """
    if not build_indices:
        pytest.skip('Skipped test as build_indices: {}'.format(build_indices))
    file_name = hisat2.HT2_FORMAT.format(index_prefix, index)
    dir_index_name = os.path.basename(os.path.normpath(dir_index))
    utils.equal_file_sizes(
        os.path.join(expected_fixture, dir_index_name, file_name),
        os.path.join(dir_index, file_name))


def compare_fq_files(expected_fixture, dir_tmp, sample, file_name):
    """
    Test FASTQ files for equality. See
    :py:func:`riboviz.fastq.equal_fastq`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: File name
    :type file_name: str or unicode
    """
    dir_tmp_name = os.path.basename(os.path.normpath(dir_tmp))
    fastq.equal_fastq(os.path.join(expected_fixture, dir_tmp_name, sample,
                                   file_name),
                      os.path.join(dir_tmp, sample,
                                   file_name))


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_cutadapt_fq(is_multiplexed, expected_fixture, dir_tmp, sample):
    """
    Test ``cutadapt`` FASTQ files for equality. See
    :py:func:`riboviz.fastq.equal_fastq`.

    Skipped if ``is_multiplexed`` is ``True``.

    :param is_multiplexed: Are the samples from multiplexed files?
    :type is_multiplexed: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if is_multiplexed:
        pytest.skip('Skipped test as is_multiplexed: {}'.format(
            is_multiplexed))
    compare_fq_files(expected_fixture, dir_tmp, sample,
                     workflow_files.ADAPTER_TRIM_FQ)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_umitools_extract_fq(extract_umis, is_multiplexed,
                             expected_fixture, dir_tmp,
                             sample):
    """
    Test ``umi_tools extract`` FASTQ files for equality. See
    :py:func:`riboviz.fastq.equal_fastq`.

    Skipped if :py:const:`riboviz.params.EXTRACT_UMIS` is ``False``
    or if ``is_multiplexed`` is ``True``.

    :param extract_umi: Configuration parameter
    :type extract_umis: bool
    :param is_multiplexed: Are the samples from multiplexed files?
    :type is_multiplexed: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if not extract_umis:
        pytest.skip('Skipped test as extract_umis: {}'.format(extract_umis))
    if is_multiplexed:
        pytest.skip('Skipped test as is_multiplexed: {}'.format(
            is_multiplexed))
    compare_fq_files(expected_fixture, dir_tmp, sample,
                     workflow_files.UMI_EXTRACT_FQ)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name", [
    workflow_files.NON_RRNA_FQ,
    workflow_files.UNALIGNED_FQ])
def test_hisat_fq(expected_fixture, dir_tmp, sample, file_name):
    """
    Test ``hisat`` FASTQ files for equality. See
    :py:func:`riboviz.fastq.equal_fastq`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_fq_files(expected_fixture, dir_tmp, sample, file_name)


def compare_sam_files(expected_directory, directory,
                      scratch_directory, sample, file_name):
    """
    Test SAM files for equality. The SAM files are sorted
    into temporary SAM files which are then compared. See
    :py:func:`riboviz.sam_bam.equal_sam`.

    :param expected_directory: Expected data directory
    :type expected_directory: str or unicode
    :param directory: Data directory
    :type directory: str or unicode
    :param scratch_directory: scratch files directory
    :type scratch_directory: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    dir_name = os.path.basename(os.path.normpath(directory))
    expected_file = os.path.join(
        expected_directory, dir_name, sample, file_name)
    actual_file = os.path.join(directory, sample, file_name)
    expected_copy_dir = os.path.join(scratch_directory, "expected")
    os.mkdir(expected_copy_dir)
    actual_copy_dir = os.path.join(scratch_directory, "actual")
    os.mkdir(actual_copy_dir)
    expected_copy_file = os.path.join(expected_copy_dir, file_name)
    actual_copy_file = os.path.join(actual_copy_dir, file_name)
    pysam.sort("-o", expected_copy_file, expected_file)
    pysam.sort("-o", actual_copy_file, actual_file)
    sam_bam.equal_sam(expected_copy_file, actual_copy_file)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name", [
    workflow_files.ORF_MAP_SAM,
    workflow_files.RRNA_MAP_SAM])
def test_hisat2_sam(expected_fixture, dir_tmp, scratch_directory,
                    sample, file_name):
    """
    Test ``hisat`` SAM files for equality. The SAM files are sorted
    into temporary SAM files which are then compared. See
    :py:func:`compare_sam_files`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param scratch_directory: scratch files directory
    :type scratch_directory: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_sam_files(expected_fixture, dir_tmp, scratch_directory,
                      sample, file_name)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_trim5p_mismatch_sam(
        trim_5p_mismatches, expected_fixture, dir_tmp,
        scratch_directory, sample):
    """
    Test :py:mod:`riboviz.tools.trim_5p_mismatch` SAM files for
    equality. The SAM files are sorted into temporary SAM files which
    are then compared. See :py:func:`compare_sam_files`.

    Skipped if :py:const:`riboviz.params.TRIM_5P_MISMATCHES` is
    ``False``.

    :param trim_5p_mismatches: Configuration parameter
    :type trim_5p_mismatches: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param scratch_directory: scratch files directory
    :type scratch_directory: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if not trim_5p_mismatches:
        pytest.skip('Skipped test as trim_5p_mismatches: {}'.format(
            trim_5p_mismatches))
    compare_sam_files(expected_fixture, dir_tmp, scratch_directory,
                      sample, workflow_files.ORF_MAP_CLEAN_SAM)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_trim5p_mismatch_tsv(
        trim_5p_mismatches, expected_fixture, dir_tmp, sample):
    """
    Test :py:mod:`riboviz.tools.trim_5p_mismatch` TSV files for
    equality. See :py:func:`riboviz.utils.equal_tsv`.

    Skipped if :py:const:`riboviz.params.TRIM_5P_MISMATCHES` is
    ``False``.

    :param trim_5p_mismatches: Configuration parameter
    :type trim_5p_mismatches: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if not trim_5p_mismatches:
        pytest.skip('Skipped test as trim_5p_mismatches: {}'.format(
            trim_5p_mismatches))
    dir_tmp_name = os.path.basename(os.path.normpath(dir_tmp))
    utils.equal_tsv(
        os.path.join(expected_fixture, dir_tmp_name, sample,
                     workflow_files.TRIM_5P_MISMATCH_TSV),
        os.path.join(dir_tmp, sample,
                     workflow_files.TRIM_5P_MISMATCH_TSV))


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_samtools_view_sort_index_orf_map_clean_bam(
        expected_fixture, dir_tmp, sample):
    """
    Test ``samtools view | samtools sort`` BAM and ``samtools index``
    BAI files for equality. See :py:func:`riboviz.sam_bam.equal_bam` and
    :py:func:`riboviz.utils.equal_file_sizes`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    dir_tmp_name = os.path.basename(os.path.normpath(dir_tmp))
    sam_bam.equal_bam(
        os.path.join(expected_fixture, dir_tmp_name, sample,
                     workflow_files.ORF_MAP_CLEAN_BAM),
        os.path.join(dir_tmp, sample,
                     workflow_files.ORF_MAP_CLEAN_BAM))
    bai_file_name = sam_bam.BAI_FORMAT.format(workflow_files.ORF_MAP_CLEAN_BAM)
    utils.equal_file_sizes(
        os.path.join(expected_fixture, dir_tmp_name, sample,
                     bai_file_name),
        os.path.join(dir_tmp, sample, bai_file_name))


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_samtools_index_dedup_bam(dedup_umis, dir_tmp, sample):
    """
    Test ``samtools index`` BAM and BAI files. Check files exist only.

    Skipped if :py:const:`riboviz.params.DEDUP_UMIS` is ``False``.

    :param dedup_umi: Configuration parameter
    :type dedup_umis: bool
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if not dedup_umis:
        pytest.skip('Skipped test as dedup_umis: {}'.format(dedup_umis))
    actual_file = os.path.join(dir_tmp, sample, workflow_files.DEDUP_BAM)
    assert os.path.exists(actual_file), "Non-existent file: %s" % actual_file
    actual_bai_file = os.path.join(
        dir_tmp, sample, sam_bam.BAI_FORMAT.format(workflow_files.DEDUP_BAM))
    assert os.path.exists(actual_bai_file),\
        "Non-existent file: %s" % actual_bai_file


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_samtools_view_sort_index(dedup_umis, expected_fixture,
                                  dir_out, sample):
    """
    Test ``samtools view | samtools sort`` BAM and ``samtools index``
    BAI files for equality. See :py:func:`riboviz.sam_bam.equal_bam`
    and :py:func:`riboviz.utils.equal_file_sizes`.

    If :py:const:`riboviz.params.DEDUP_UMIS` is ``False`` then
    only the existence of the files are checked as these files can
    differ between runs depending on which reads are removed by
    ``umi_tools dedup``.

    :param dedup_umi: Configuration parameter
    :type dedup_umis: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    file_name = sam_bam.BAM_FORMAT.format(sample)
    bai_file_name = sam_bam.BAI_FORMAT.format(file_name)
    dir_out_name = os.path.basename(os.path.normpath(dir_out))
    expected_file = os.path.join(
        expected_fixture, dir_out_name, sample, file_name)
    actual_file = os.path.join(dir_out, sample, file_name)
    assert os.path.exists(actual_file), "Non-existent file: %s" % actual_file
    expected_bai_file = os.path.join(
        expected_fixture, dir_out_name, sample, bai_file_name)
    actual_bai_file = os.path.join(dir_out, sample, bai_file_name)
    assert os.path.exists(actual_bai_file),\
        "Non-existent file: %s" % actual_bai_file
    if dedup_umis:
        return
    sam_bam.equal_bam(expected_file, actual_file)
    utils.equal_file_sizes(expected_bai_file, actual_bai_file)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("stats_file", ["edit_distance.tsv",
                                        "per_umi_per_position.tsv",
                                        "per_umi.tsv"])
def test_umitools_dedup_stats_tsv(
        dedup_umis, dedup_stats, dir_tmp, sample, stats_file):
    """
    Test ``umi_tools dedup --output-stats`` TSV files exist.

    As these files can differ between runs depending on which reads
    are removed by ``umi_tools dedup``, only the existence of the
    files is checked.

    Skipped if :py:const:`riboviz.params.DEDUP_UMIS` is ``False``
    or :py:const:`riboviz.params.DEDUP_STATS` is ``False``.

    :param dedup_umi: Configuration parameter
    :type dedup_umis: bool
    :param dedup_stats: Configuration parameter
    :type dedup_stats: bool
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param stats_file: statistics file name
    :type stats_file: str or unicode
    """
    if not dedup_umis:
        pytest.skip('Skipped test as dedup_umis: {}'.format(dedup_umis))
    if not dedup_stats:
        pytest.skip('Skipped test as dedup_stats: {}'.format(dedup_stats))
    file_name = os.path.join(sample,
                             workflow_files.DEDUP_STATS_FORMAT.format(
                                 stats_file))
    actual_file = os.path.join(dir_tmp, file_name)
    assert os.path.exists(actual_file), "Non-existent file: %s" % actual_file


def compare_tsv_files(expected_fixture, directory, sample, file_name):
    """
    Test TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param directory: Actual data directory
    :type directory: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    directory_name = os.path.basename(os.path.normpath(directory))
    expected_file = os.path.join(expected_fixture, directory_name,
                                 sample, file_name)
    utils.equal_tsv(
        expected_file,
        os.path.join(directory, sample, file_name))


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_umitools_pre_dedup_group_tsv(
        dedup_umis, group_umis, expected_fixture, dir_tmp, sample):
    """
    Test ``umi_tools group`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    Skipped if :py:const:`riboviz.params.DEDUP_UMIS` is ``False``
    or :py:const:`riboviz.params.GROUP_UMIS` is ``False``.

    :param dedup_umis: Configuration parameter
    :type dedup_umis: bool
    :param group_umis: Configuration parameter
    :type group_umis: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if not dedup_umis:
        pytest.skip('Skipped test as dedup_umis: {}'.format(dedup_umis))
    if not group_umis:
        pytest.skip('Skipped test as group_umis: {}'.format(group_umis))
    compare_tsv_files(expected_fixture, dir_tmp, sample,
                      workflow_files.PRE_DEDUP_GROUPS_TSV)


@pytest.mark.usefixtures("skip_index_tmp_fixture")
@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_umitools_post_dedup_group_tsv(
        dedup_umis, group_umis, dir_tmp, sample):
    """
    Test ``umi_tools group`` TSV file exists.

    As these files can differ between runs depending on which reads
    are removed by ``umi_tools dedup``, only the existence of the file
    is checked.

    Skipped if :py:const:`riboviz.params.DEDUP_UMIS` is ``False``
    or :py:const:`riboviz.params.GROUP_UMIS` is ``False``.

    :param dedup_umis: Configuration parameter
    :type dedup_umis: bool
    :param group_umis: Configuration parameter
    :type group_umis: bool
    :param dir_tmp: Temporary directory
    :type dir_tmp: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    if not dedup_umis:
        pytest.skip('Skipped test as dedup_umis: {}'.format(dedup_umis))
    if not group_umis:
        pytest.skip('Skipped test as group_umis: {}'.format(group_umis))
    actual_file = os.path.join(
        dir_tmp, sample, workflow_files.POST_DEDUP_GROUPS_TSV)
    assert os.path.exists(actual_file), "Non-existent file: %s" % actual_file


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name", [
    workflow_files.MINUS_BEDGRAPH,
    workflow_files.PLUS_BEDGRAPH])
def test_bedtools_bedgraph(expected_fixture, make_bedgraph, dir_out,
                           sample, file_name):
    """
    Test ``bedtools genomecov`` bedgraph files for equality. See
    :py:func:`riboviz.bedgraph.equal_bedgraph`.

    Skipped if :py:const:`riboviz.params.MAKE_BEDGRAPH` is ``False``.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param make_bedgraph: Configuration parameter
    :type make_bedgraph: bool
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not make_bedgraph:
        pytest.skip('Skipped test as make_bedgraph: {}'.format(make_bedgraph))
    dir_out_name = os.path.basename(os.path.normpath(dir_out))
    expected_file = os.path.join(expected_fixture, dir_out_name,
                                 sample, file_name)
    bedgraph.equal_bedgraph(expected_file,
                            os.path.join(dir_out, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_bam_to_h5_h5(expected_fixture, dir_out, sample):
    """
    Test ``bam_to_h5.R`` H5 files for equality. See
    :py:func:`riboviz.h5.equal_h5`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    """
    file_name = h5.H5_FORMAT.format(sample)
    dir_out_name = os.path.basename(os.path.normpath(dir_out))
    expected_file = os.path.join(expected_fixture, dir_out_name,
                                 sample, file_name)
    h5.equal_h5(expected_file,
                os.path.join(dir_out, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.ORF_TPMS_AND_COUNTS_TSV,
                          workflow_r.METAGENE_START_STOP_READ_COUNTS_TSV,
                          workflow_r.METAGENE_NORMALIZED_PROFILE_START_STOP_TSV,
                          workflow_r.READ_COUNTS_BY_LENGTH_TSV,
                          workflow_r.METAGENE_POSITION_LENGTH_COUNTS_TSV])
def test_generate_stats_figs_tsv(expected_fixture, dir_out, sample,
                                 file_name):
    """
    Test ``generate_stats_figs.R`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    compare_tsv_files(expected_fixture, dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.NT_FREQ_PER_READ_POSITION_TSV])
def test_generate_stats_figs_metagene_tsv(
        output_metagene_normalized_profile,
        expected_fixture, dir_out, sample, file_name):
    """
    Test ``generate_stats_figs.R`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    Skipped if
    :py:const:`riboviz.params.OUTPUT_METAGENE_NORMALIZED_PROFILE` is
    ``False``.

    :param output_metagene_normalized_profile: Configuration parameter
    :type output_metagene_normalized_profile: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not output_metagene_normalized_profile:
        pytest.skip(
            'Skipped test as output_metagene_normalized_profile: {}'.format(
                output_metagene_normalized_profile))
    compare_tsv_files(expected_fixture, dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.ORF_TPMS_VS_FEATURES_TSV])
def test_generate_stats_figs_features_tsv(
        features_file, expected_fixture, dir_out, sample, file_name):
    """
    Test ``generate_stats_figs.R`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    Skipped if :py:const:`riboviz.params.FEATURES_FILE` is
    ``None``.

    :param features_file: Configuration parameter
    :type features_file: str or unicode
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not features_file:
        pytest.skip('Skipped test as features_file: {}'.format(
            features_file))
    compare_tsv_files(expected_fixture, dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.NORMALIZED_DENSITY_APESITES_PER_CODON_TSV,
                          workflow_r.NORMALIZED_DENSITY_APESITES_PER_CODON_LONG_TSV])
def test_generate_stats_figs_t_rna_codon_positions_tsv(
        t_rna_file, codon_positions_file, expected_fixture, dir_out,
        sample, file_name):
    """
    Test ``generate_stats_figs.R`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    Skipped if :py:const:`riboviz.params.T_RNA_FILE` is ``None`` or
    :py:const:`riboviz.params.CODON_POSITIONS_FILE` is ``None``.

    :param t_rna_file: Configuration parameter
    :type t_rna_file: bool
    :param codon_positions_file: Configuration parameter
    :type codon_positions_file: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not t_rna_file:
        pytest.skip('Skipped test as t_rna_file: {}'.format(
            t_rna_file))
    if not codon_positions_file:
        pytest.skip('Skipped test as codon_positions_file: {}'.format(
            codon_positions_file))
    compare_tsv_files(expected_fixture, dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.READ_FRAME_PER_ORF_TSV,
                          workflow_r.READ_FRAME_PER_ORF_FILTERED_TSV])
def test_generate_stats_figs_asite_disp_length_tsv(
        asite_disp_length_file, expected_fixture, dir_out, sample,
        file_name):
    """
    Test ``generate_stats_figs.R`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    Skipped if :py:const:`riboviz.params.ASITE_DISP_LENGTH_FILE` is
    ``None``.

    :param asite_disp_length_file: Configuration parameter
    :type asite_disp_length_file: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not asite_disp_length_file:
        pytest.skip('Skipped test as asite_disp_length_file: {}'.format(
            asite_disp_length_file))
    compare_tsv_files(expected_fixture, dir_out, sample, file_name)


def check_pdf_file_exists(dir_out, sample, file_name):
    """
    Check that a PDF file exists.

    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    :raise AssertionError: if PDF file does not exist
    """
    actual_file = os.path.join(dir_out, sample, file_name)
    assert os.path.exists(actual_file), "Non-existent file: %s" % actual_file


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.READ_COUNTS_BY_LENGTH_PDF,
                          workflow_r.METAGENE_START_STOP_READ_COUNTS_PDF,
                          workflow_r.METAGENE_START_BARPLOT_BY_LENGTH_PDF,
                          workflow_r.METAGENE_START_RIBOGRID_BY_LENGTH_PDF,
                          workflow_r.METAGENE_NORMALIZED_PROFILE_START_STOP_PDF])
def test_generate_stats_figs_pdf(
        output_pdfs, dir_out, sample, file_name):
    """
    Test ``generate_stats_figs.R`` PDF files exist.

    Skipped if :py:const:`riboviz.params.OUTPUT_PDFS` is ``False``.

    :param output_pdfs: Configuration parameter
    :type output_pdfs: bool
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not output_pdfs:
        pytest.skip('Skipped test as output_pdfs: {}'.format(output_pdfs))
    check_pdf_file_exists(dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.ORF_TPMS_VS_FEATURES_PDF])
def test_generate_stats_figs_features_pdf(
        features_file, output_pdfs, dir_out, sample, file_name):
    """
    Test ``generate_stats_figs.R`` PDF files exist.

    Skipped if :py:const:`riboviz.params.FEATURES_FILE` is
    ``None``.

    Skipped if :py:const:`riboviz.params.OUTPUT_PDFS` is ``False``.

    :param features_file: Configuration parameter
    :type features_file: book
    :param output_pdfs: Configuration parameter
    :type output_pdfs: bool
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not features_file:
        pytest.skip('Skipped test as features_file: {}'.format(
            features_file))
    if not output_pdfs:
        pytest.skip('Skipped test as output_pdfs: {}'.format(output_pdfs))
    check_pdf_file_exists(dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.NORMALIZED_DENSITY_APESITES_PER_CODON_PDF])
def test_generate_stats_figs_t_rna_codon_positions_pdf(
        t_rna_file, codon_positions_file, output_pdfs, dir_out,
        sample, file_name):
    """
    Test ``generate_stats_figs.R`` PDF files exist.

    Skipped if :py:const:`riboviz.params.T_RNA_FILE` is ``None`` or
    :py:const:`riboviz.params.CODON_POSITIONS_FILE` is ``None``.

    Skipped if :py:const:`riboviz.params.OUTPUT_PDFS` is ``False``.

    :param t_rna_file: Configuration parameter
    :type t_rna_file: bool
    :param codon_positions_file: Configuration parameter
    :type codon_positions_file: bool
    :param output_pdfs: Configuration parameter
    :type output_pdfs: bool
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not t_rna_file:
        pytest.skip('Skipped test as t_rna_file: {}'.format(
            t_rna_file))
    if not codon_positions_file:
        pytest.skip('Skipped test as codon_positions_file: {}'.format(
            codon_positions_file))
    if not output_pdfs:
        pytest.skip('Skipped test as output_pdfs: {}'.format(output_pdfs))
    check_pdf_file_exists(dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name",
                         [workflow_r.FRAME_PROPORTIONS_PER_ORF_PDF])
def test_generate_stats_figs_asite_disp_length_pdf(
        asite_disp_length_file, output_pdfs, dir_out, sample,
        file_name):
    """
    Test ``generate_stats_figs.R`` PDF files exist.

    Skipped if :py:const:`riboviz.params.ASITE_DISP_LENGTH_FILE` is
    ``None``.

    Skipped if :py:const:`riboviz.params.OUTPUT_PDFS` is ``False``.

    :param asite_disp_length_file: Configuration parameter
    :type asite_disp_length_file: bool
    :param output_pdfs: Configuration parameter
    :type output_pdfs: bool
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not asite_disp_length_file:
        pytest.skip('Skipped test as asite_disp_length_file: {}'.format(
            asite_disp_length_file))
    if not output_pdfs:
        pytest.skip('Skipped test as output_pdfs: {}'.format(output_pdfs))
    check_pdf_file_exists(dir_out, sample, file_name)


@pytest.mark.usefixtures("prep_riboviz_fixture")
@pytest.mark.parametrize("file_name", [workflow_files.STATIC_HTML_FILE])
def test_analysis_outputs_html(run_static_html, expected_fixture,
                               dir_out, sample, file_name):
    """
    Test ``AnalysisOutputs.Rmd`` html files for equality. See
    :py:func:`riboviz.html.equal_html`.

    Skipped if :py:const:`riboviz.params.RUN_STATIC_HTML` is ``False``.

    :param run_static_html: Configuration parameter
    :type run_static_html: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    :param sample: Sample name
    :type sample: str or unicode
    :param file_name: file name
    :type file_name: str or unicode
    """
    if not run_static_html:
        pytest.skip('Skipped test as run_static_html: {}'.format(
            run_static_html))
    file_name = workflow_files.STATIC_HTML_FILE.format(sample)
    dir_out_name = os.path.basename(os.path.normpath(dir_out))
    expected_file = os.path.join(expected_fixture, dir_out_name,
                                 sample, file_name)
    assert os.path.exists(os.path.join(dir_out, sample, file_name))
    html.equal_html(expected_file,
                    os.path.join(dir_out, sample, file_name))


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_collate_orf_tpms_and_counts_tsv(expected_fixture, dir_out):
    """
    Test ``collate_tpms.R`` TSV files for equality. See
    :py:func:`riboviz.utils.equal_tsv`.

    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    """
    dir_out_name = os.path.basename(os.path.normpath(dir_out))
    expected_file = os.path.join(expected_fixture, dir_out_name,
                                 workflow_r.TPMS_ALL_CDS_ALL_SAMPLES_TSV)
    utils.equal_tsv(expected_file,
                    os.path.join(dir_out, workflow_r.TPMS_ALL_CDS_ALL_SAMPLES_TSV),
                    ignore_row_order=True,
                    na_to_empty_str=True)


@pytest.mark.usefixtures("prep_riboviz_fixture")
def test_read_counts_per_file_tsv(count_reads, expected_fixture, dir_out):
    """
    Test :py:mod:`riboviz.tools.count_reads` TSV files for
    equality. See :py:func:`riboviz.count_reads.equal_read_counts`.

    Skipped if :py:const:`riboviz.params.COUNT_READS` is ``False``.

    :param count_reads: Configuration parameter
    :type count_reads: bool
    :param expected_fixture: Expected data directory
    :type expected_fixture: str or unicode
    :param dir_out: Output directory
    :type dir_out: str or unicode
    """
    if not count_reads:
        pytest.skip('Skipped test as count_reads: {}'.format(count_reads))
    dir_out_name = os.path.basename(os.path.normpath(dir_out))
    count_reads_module.equal_read_counts(
        os.path.join(expected_fixture, dir_out_name,
                     workflow_files.READ_COUNTS_PER_FILE_FILE),
        os.path.join(dir_out, workflow_files.READ_COUNTS_PER_FILE_FILE))
