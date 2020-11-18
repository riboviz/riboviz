# Memory and storage

Information and advice relating to RiboViz's memory and storage requirements.

---

## Memory requirements

The memory requirements of RiboViz depend upon the datasets being processed.

The RiboViz developers use development environments with 8-16GB RAM. This has been found to be adequate to handle the vignette - [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md). Multiple full size yeast datasets have been successfully processed on developers' own machines with 16GB RAM available.

One developer had issues when running the vignette as batch job on the University of Edinburgh ECDF Linux Compute Cluster) [Eddie](https://www.ed.ac.uk/information-services/research-support/research-computing/ecdf/high-performance-computing). Their job requested an 8GB node but this was terminated as the job exceeded this. On requesting a 16GB node, the job ran to success, with the system reporting that almost 14GB has been used. It is unclear at present as to why this is.

Full size yeast datasets have been successfully processed on EDDIE with 16GB RAM.

### Troubleshooting: deduplication and memory issues

If running a workflow using deduplication, then you might encounter memory issues arising from an issue with UMI-Tools. UMI-tools has an option, `--output-stats`, which calculates statistics relating to deduplication. This option can result in excessive memory usage during deduplication (for more information, see the UMI-Tools issue "excessive dedup memory usage with output-stats" [#409](https://github.com/CGATOxford/UMI-tools/issues/409)).

The workflow sets `--output-stats` by default. If you find you are having issues with memory usage then you might be able to resolve these by requesting that deduplication statistics are not calculated. You can do this by editing your YAML configuration file and setting:

```yaml
dedup_stats: FALSE
```

### Profiling RiboViz's memory requirements

We have a ticket for a future release to Profile RiboViz to determine memory requirements [#179](https://github.com/riboviz/riboviz/issues/179). In the meantime, welcome contributions from users as to the memory requirements of your analyses, both when running the vignette and your own analyses. Please record an indication of the hardware and operating system you used, the memory you had available and the size of your input files (FASTA, GFF ans FASTQ files) (both bytes and number of lines).

---

## Storage requirements

The workflow generates many intermediate files and some of these may be unompressed and **large**, i.e. about the same size as the input files. All these files are placed in a temporary directory (`dir_tmp`). The temporary directory's contents can be inspected for troubleshooting, if necessary.

The Python workflow also creates numerous log files.

The Nextflow `work` directory has the originals of all temporary and output files.

For example, here is the volume of the outputs from a run of the vignette as documented in [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md):

* Python workflow:

| Directory         |   MB |
| ----------------- | ---- |
| `vignette/index`  |    9 |
| `vignette/tmp`    |  812 |
| `vignette/output` |    3 |
| `vignette/logs`   |    1 |
| Total             |  825 |

* Nextflow workflow:

| Directory         |   MB |
| ----------------- | ---- |
| `vignette/index`  |    9 |
| `vignette/tmp`    |  813 |
| `vignette/output` |    3 |
| `work`            |  826 |
| Total             | 1651 |

**Tip:** We recommend you regularly delete temporary directories, log directories (Python workflow) and work directory (Nextflow workflow) when you have completed your analyses to your satisfaction.
