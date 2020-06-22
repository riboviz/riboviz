# Memory and storage usage

## Managing your storage usage

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
