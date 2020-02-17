"""
Helper methods for comparing and validating files of different types.
"""
import os
import os.path
from riboviz import bedgraph
from riboviz import fastq
from riboviz import h5
from riboviz import hisat2
from riboviz import sam_bam
from riboviz import utils


def compare_files(file1, file2, compare_names=True):
    """
    Compare two files for equality. The following functions are used
    to compare each type of file:

    * pdf: utils.equal_file_ames(file1, file2)
    * ht2, .bai: utils.equal_file_sizes(file1, file2)
    * h5: h5.equal_h5(file1, file2)
    * bedgraph: bedgraph.equal_bedgraph(file1, file2)
    * bam: sam_bam.equal_bam(file1, file2)
    * sam: sam_bam.equal_sam(file1, file2)
    * tsv: utils.equal_tsv(file1, file2)
    * fq: fastq.equal_fastq(file1, file2)

    :param file1: File name
    :type file1: str or unicode
    :param file2: File name
    :type file2: str or unicode
    :param compare_names: compare file names?
    :type: bool
    :raise AssertionError: if files differ
    """
    assert os.path.exists(file1), "Non-existent file: %s" % file1
    assert os.path.exists(file2), "Non-existent file: %s" % file2
    assert not os.path.isdir(file1), "Directory: %s" % file1
    assert not os.path.isdir(file2), "Directory: %s" % file2
    if compare_names:
        utils.equal_file_names(file1, file2)
    ext = utils.get_file_ext(file1)
    if ext in ["pdf"]:
        utils.equal_file_names(file1, file2)
    elif ext in [hisat2.HT2_EXT, sam_bam.BAI_EXT]:
        utils.equal_file_sizes(file1, file2)
    elif ext in [h5.H5_EXT]:
        h5.equal_h5(file1, file2)
    elif ext in [bedgraph.BEDGRAPH_EXT]:
        bedgraph.equal_bedgraph(file1, file2)
    elif ext in [sam_bam.BAM_EXT]:
        sam_bam.equal_bam(file1, file2)
    elif ext in [sam_bam.SAM_EXT]:
        sam_bam.equal_sam(file1, file2)
    elif ext in ["tsv"]:
        utils.equal_tsv_files(file1, file2)
    elif ext in fastq.FASTQ_ALL_EXTS:
        fastq.equal_fastq(file1, file2)
