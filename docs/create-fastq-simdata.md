# create_fastq_simdata.py example FASTQ file generator

`riboviz/tools/create_fastq_simdata.py` creates simple simulated FASTQ files to test adaptor trimming, UMI extraction and deduplication. The files in `data/example/` were created using this script.

---

## Usage

To run:

```console
$ python -m riboviz.tools.create_fastq_simdata <DIRECTORY>
```

Or:

```console
$ python -m riboviz/tools/create_fastq_simdata.py <DIRECTORY>
```

where `<DIRECTORY>` is the directory into which the simulated files are to be written.

---

## Example simulated data files

The following files are created.

### Reads with 5' UMI, 3' UMI and adaptor

These files can be used to test adaptor trimming and deduplication.

```
umi5_umi3_umi_adaptor.fastq
```

FASTQ file with 9 reads, each with a 4nt UMI at the 5' end, a 4nt UMI at the 3' end and a 11nt adaptor at the 3' end. Reads can be grouped by UMI into 5 groups.

```
umi5_umi3_umi.fastq
````

FASTQ file identical to the above but with the adaptor trimmed.

```
umi5_umi3.fastq
```

FASTQ file identical to the above but with the UMIs extracted and concatenated to the header, with a "_" delimiter.

### Reads with 3' UMI and adaptor

These files can be used to test adaptor trimming and deduplication.

```
umi3_umi_adaptor.fastq
```

FASTQ file with 8 reads, each with a 4nt UMI at the 3' end and a 11nt adaptor at the 3' end. Reads can be grouped by UMI into 4 groups.

```
umi3_umi.fastq
```

FASTQ file identical to the above but with the adaptor trimmed.

```
umi3.fastq`
```

FASTQ file identical to the above but with the UMI extracted and concatenated to the header, with a "_" delimiter.

### Reads with 5' UMI, 3' UMI, adaptor and barcode

These files can be used to test adaptor trimming, demultiplexing and deduplication.

```
multiplex_barcodes.tsv
```

Tab-separated values file with `SampleID` column (with values `Tag0|1|2`) and `TagRead` column (with values `ACG`, `GAC`, `CGA`). This is consistent with the file format expected by `riboviz.tools.demultiplex_fastq`.

```
multiplex_umi_barcode_adaptor.fastq
```

FASTQ file with 90 reads:

* Each read has a 4nt UMI at the 5' end, a 4nt UMI at the 3' end, a 3nt barcode at the 3' end and a 11nt adaptor at the 3' end.
* There are 9 reads for each of the following barcodes:
  - `ACG`, `ACT`, `TAG`
  - `GAC`, `GTC`, `GTA`
  - `CGA`, `TGA`, `CTT`
* The second and third barcodes in each list have a mismatch of 1nt and 2nt respectively with the first barcode in each list.
* When the file is demultiplexed, assuming up to 2 mismatches are allowed, then 3 sets of 27 reads will be produced, grouped by the 1st barcode in each list.
* There are 9 reads with barcode `TTT`, which has a mismatch of 3nts to `ACG`, `GAC`, `CGA`. When the file is demultiplexed, assuming up to 2 mismatches are allowed, then these 9 reads will be unassigned.

```
multiplex_umi_barcode.fastq
```

FASTQ file identical to the above but with the adaptor trimmed.

```
multiplex.fastq
```

FASTQ file identical to the above but with the barcode and UMIs extracted into the header and delimited by "_".

```
deplex/multiplex_tag0|1|2.fastq
```

FASTQ files each with 27 reads representing the results expected when demultiplexing `multiplex.fastq` using `riboviz.tools.demultiplex_fastq` and `multiplex_barcodes.tsv`.

```
deplex/Unassigned.fastq
```

FASTQ files with 9 reads representing the unassigned reads (those with barcode `TTT`) expected when demultiplexing `multiplex.fastq` using `riboviz.tools.demultiplex_fastq` and `multiplex_barcodes.tsv`.

```
deplex/num_reads.tsv
```

Tab-separated values with expected counts of reads for each barcode expected when demultiplexing `multiplex.fastq` using `riboviz.tools.demultiplex_fastq` and `multiplex_barcodes.tsv`.
