# demultiplex-fastq.py fastq demultiplexer

`riboviz/tools/demultiplex-fastq.py` is a command-line tool to demultiplex fastq files using UMI-tools-compliant barcodes present within the fastq headers. These headers are assumed to be of form:

```
@..._<BARCODE>_...
```

where the barcode is the first section delimited by underscores. If another delimiter was used then that can be specified.

---

## Usage

```
$ python -m riboviz.tools.demultiplex_fastq  -h
usage: demultiplex_fastq.py [-h] [-ss [SAMPLE_SHEET_FILE]] [-r1 [READ1_FILE]]
                            [-r2 [READ2_FILE]] [-m MISMATCHES] [-o [OUT_DIR]]
                            [-d [DELIMITER]]

Demultiplex reads from fastq[.gz] by inline barcodes

optional arguments:
  -h, --help            show this help message and exit
  -ss [SAMPLE_SHEET_FILE], --samplesheet [SAMPLE_SHEET_FILE]
                        Sample sheet filename, tab-delimited text format with
                        SampleID and TagRead columns
  -r1 [READ1_FILE], --read1 [READ1_FILE]
                        Read 1 filename, fastq[.gz] format
  -r2 [READ2_FILE], --read2 [READ2_FILE]
                        Read 2 pair filename, fastq[.gz] format
  -m MISMATCHES, --mismatches MISMATCHES
                        Number of mismatches permitted in barcode
  -o [OUT_DIR], --outdir [OUT_DIR]
                        Output directory
  -d [DELIMITER], --delimiter [DELIMITER]
                        Barcode delimiter
```

Inputs:

* `-ss|--samplesheet`: Sample sheet filename, tab-delimited text format with SampleID and TagRead (barcode) columns
* `-r1|--read1`: Read 1 filename, fastq[.gz] format
* `-r2|--read2`: Read 2 pair filename, fastq[.gz] format (must be consistent with Read 1 filename) (optional) . If provided then the read files should have read pairs in corresponding positions.
* `-m|--mismatches`: Number of mismatches permitted in barcode (optional, default 1)
* `-o|--outdir`: Output directory (optional, default output)
* `-d|--delimiter`: Barcode delimiter (optional, default "_")

Outputs:

* If `r1` only was provided:
  - A file `SampleID.fastq[.gz]` with assigned reads.
  - A file, `Unassigned.fastq[.gz]`, with information on unassigned reads.
* If `r1` and `r2` were provided:
  - Files `SampleID_R1.fastq[.gz]` and `SampleID_R2.fastq[.gz]` with assigned reads.
  - Files, `Unassigned_R1.fastq[.gz]` and `Unassigned_R2.fastq[.gz]` with information on unassigned reads.
* If the input file(s) had were of type fastq.gz then the output files will be of type fastq.gz.
* A file, `num_reads.tsv`, with `SampleID`, `TagRead` and `NumReads` columns, specifying the number of reads for each `SampleID` and `TagRead` in the original sample sheet, plus the number of unassigned reads and the total number of reads.

---

## Example:

Run UMI-tools on sample data and extract barcodes:

```console
$ mkdir extracts/
$ umi_tools extract --bc-pattern="^(?P<cell_1>.{9})(?P<umi_1>.{0}).+$" \
  --extract-method=regex -I data/demultiplex/Sample_4reads_R1.fastq.gz \
  -S extracts/Sample_4reads_R1.fastq.gz
$ umi_tools extract --bc-pattern="^(?P<cell_1>.{9})(?P<umi_1>.{0}).+$" \
  --extract-method=regex -I data/demultiplex/Sample_4reads_R2.fastq.gz \
  -S extracts/Sample_4reads_R2.fastq.gz
$ umi_tools extract --bc-pattern="^(?P<cell_1>.{9})(?P<umi_1>.{0}).+$" \
  --extract-method=regex -I data/demultiplex/Sample_init10000_R1.fastq.gz \
  -S extracts/Sample_init10000_R1.fastq.gz
$ umi_tools extract --bc-pattern="^(?P<cell_1>.{9})(?P<umi_1>.{0}).+$" \
  --extract-method=regex -I data/demultiplex/Sample_init10000_R2.fastq.gz \
  -S extracts/Sample_init10000_R2.fastq.gz
```

(note that the regular expression to extract the barcodes intentionally includes a zero-length UMI - we are assuming that this data has no UMIs)

Demultiplex single-end data:

```console
$ mkdir extracts-deplexed/
$ python -m riboviz.tools.demultiplex_fastq \
  -r1 extracts/Sample_4reads_R1.fastq.gz \
  -ss data/demultiplex/TagSeqBarcodedOligos2015.txt \
  -o extracts-deplexed/TestSingleSplit4reads
$ python -m riboviz.tools.demultiplex_fastq \
  -r1 extracts/Sample_init10000_R1.fastq.gz \
  -ss data/demultiplex/TagSeqBarcodedOligos2015.txt \
  -o extracts-deplexed/TestSingleSplit10000
```

Demultiplex paired-end data:

```console
$ python -m riboviz.tools.demultiplex_fastq \
  -r1 extracts/Sample_4reads_R1.fastq.gz \
  -r2 extracts/Sample_4reads_R2.fastq.gz \
  -ss data/demultiplex/TagSeqBarcodedOligos2015.txt \
  -o extracts-deplexed/TestPairSplit4reads
$ python -m riboviz.tools.demultiplex_fastq \
  -r1 extracts/Sample_init10000_R1.fastq.gz \
  -r2 extracts/Sample_init10000_R2.fastq.gz \
  -ss data/demultiplex/TagSeqBarcodedOligos2015.txt \
  -o extracts-deplexed/TestPairSplit10000
```

---

## Known issue - mismatches and barcode Hamming distances

If the number of mismatches (`-m|--mismatches`) is less than the Hamming distance between the barcodes (`TagReads` within the sample sheet) then a read will be assigned to the first barcode that matches even if this is not the closest barcode in terms of Hamming distance.

For example, imagine we had a barcode in a read, AGA, and our barcodes in our samplesheet are AAA, CCC, GGG, TTT. The Hamming distances between the barcode and the sample barcodes are as follows:

* d(AGA, AAA) = 1
* d(AGA, GGG) = 2
* d(AGA, TTT) = 3
* d(AGA, CCC) = 3

If mismatches is 2 then AGA could be assigned to AAA or GGG, depending on the ordering of barcodes in the sample sheet, even though AAA is closest in terms of Hamming distance.

If mismatches is 3 then AGA could be assigned to AAA, GGG, TTT depending on the ordering of barcodes in the sample sheet.

Caution should be taken if the Hamming distance of the barcodes in the sample sheet is less than the number of mismatches times 2. In the above two examples, the Hamming distance between each of the sample barcodes is 3 is less than the number of mismatches times 2, which is 4 and 6 respectively.
