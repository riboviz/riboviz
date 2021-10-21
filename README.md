# riboviz

**Ribosome profiling** provides a detailed global snapshot of protein synthesis in a cell.  At its core, this technique makes use of the observation that a translating ribosome protects around 30 nucleotides of the mRNA from nuclease activity.  High-throughput sequencing of these ribosome protected fragments (called ribosome footprints) offers a precise record of the number and location of the ribosomes at the time at which translation is stopped. Mapping the position of the ribosome protected fragments indicates the translated regions within the transcriptome.  Moreover, ribosomes spend different periods of time at different positions, leading to variation in the footprint density along the mRNA transcript. This provides an estimate of how much protein is being produced from each mRNA. Importantly, ribosome profiling is as precise and detailed as RNA sequencing. Even in a short time, since its introduction in 2009, ribosome profiling has been playing a key role in driving biological discovery.

We have developed this bioinformatics toolkit, **riboviz**, for analyzing ribosome profiling datasets. **riboviz** consists of a comprehensive and flexible analysis pipeline. The current version of **riboviz** is designed for yeast datasets.

Existing yeast datasets consist of a mix of studies, some of which use elongation inhibitors such as cycloheximide (CHX) and others that flash freeze (FF) the samples to prevent initiation and elongation during sample preparation. In general, each experimental step can potentially introduce biases in processed datasets. **riboviz** can help identify these biases by allowing users to compare and contrast datasets obtained under different experimental conditions.

The current version of **riboviz** is designed for yeast datasets but it can be customised and used to analyse datasets relating to other organisms.

For information on **riboviz**, see "riboviz: analysis and visualization of ribosome profiling datasets", Carja et al., BMC Bioinformatics 2017. doi:10.1186/s12859-017-1873-8.

All the code for processing the raw reads is available in this repository.

## Use riboviz

Quick start:

* [Install riboviz and dependencies](./docs/user/install.md)
* [Quick install scripts](./docs/user/quick-install.md) (Ubuntu and CentOS only)
* [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./docs/user/run-vignette.md). Run a "vignette" of the **riboviz** workflow to see **riboviz**'s capabilities.
* [Run an example of a configuration using environment variable tokens](./docs/user/run-env-vars-example.md).
* [Run UMI extraction, deduplication and demultiplexing examples](./docs/user/run-dedup-demultiplex-examples.md). Run **riboviz** on simulated data, to see how **riboviz** handles duplicated and multiplexed data.
* [Upgrade configuration files to current version](./docs/user/upgrade-config.md)

Usage:

* [What the riboviz workflow does](./docs/user/prep-riboviz-operation.md)
* [Configuring the riboviz workflow](./docs/user/prep-riboviz-config.md)
* [Generate YAML configuration file](https://riboviz.shinyapps.io/generate-yaml/). Use an online tool to generate a riboviz YAML configuration file.
* [Running the riboviz Nextflow workflow](./docs/user/prep-riboviz-run-nextflow.md)  
* [How To Run the riboviz Interactive Data Visualization On Your Data](docs/user/run-run_shiny_server-operation.md)
* [riboviz output files and figures](./docs/user/riboviz-outputs.md)
* [Running the riboviz workflow on Eddie](./docs/user/run-on-eddie.md)
* [Memory and storage](./docs/user/memory-storage.md). Information and advice relating to **riboviz**'s memory and storage requirements.

## Develop riboviz

General:

* [Git branching model](./docs/developer/git-branching-model.md)
* [Style guide](./docs/developer/style-guide.md): Style guidelines for the **riboviz** source code, documentation, parameters and files.

Development:

* [Install developer dependencies](./docs/developer/install.md)
* [Developing Python components](./docs/developer/dev-python.md)
* [Developing R components](./docs/developer/dev-r.md)
* [Developing Nextflow workflow](./docs/developer/dev-nextflow.md)
* [Adding, using, renaming, and removing configuration parameters](./docs/developer/config.md)
* [Adding, renaming, and removing temporary or output files](./docs/developer/output-files.md)
* [Adding and updating dependencies](./docs/developer/dependencies.md)
* [Developing and running integration tests](./docs/developer/integration-tests.md)
* [Writing and updating documentation](./docs/developer/documentation.md)

Releasing:

* [Creating a test data repository](./docs/developer/create-test-data-repository.md)
* [Building a release](./docs/developer/releasing.md)

## Reference

* [riboviz command-line tools](docs/user/command-line-tools.md).
* [Content and provenance of repository data files](./docs/reference/data.md).
* [Structure of HDF5 data](./docs/reference/hdf5-data.md). Structure of HDF5 data files produced by `bam_to_h5.R`.

## Releases

| Release | Description |
| ------- | ----------- |
| [2.1](https://github.com/riboviz/riboviz/releases/tag/2.1) | 2.1, current stable release |
| [2.0](https://github.com/riboviz/riboviz/releases/tag/2.0) | 2.0 |
| [2.0.beta](https://github.com/riboviz/riboviz/releases/tag/2.0.beta) | 2.0 beta release |
| [1.1.0](https://github.com/riboviz/riboviz/releases/tag/1.1.0) | Most recent version prior to commencement of BBSRC/NSF riboviz project |
| [1.0.0](https://github.com/riboviz/riboviz/releases/tag/1.0.0) | Associated with Carja et al. (2017) "riboviz: analysis and visualization of ribosome profiling datasets", BMC Bioinformatics, volume 18, article 461 (2017), 25 October 2017, doi: [10.1186/s12859-017-1873-8](https://doi.org/10.1186/s12859-017-1873-8) |
| [0.9.0](https://github.com/riboviz/riboviz/releases/tag/0.9.0) | Additional code/data associated with the paper below |
| [0.8.0](https://github.com/riboviz/riboviz/releases/tag/0.8.0) | Associated with Carja et al. (2017) "riboviz: analysis and visualization of ribosome profiling datasets", bioRXiv, 12 January 2017,doi: [10.1101/100032](https://doi.org/10.1101/100032) |

## Citing riboviz

To cite **riboviz**, please use both of the following references:

Cope AL, Anderson F, Favate J, Jackson M, Mok A, Kurowska A, MacKenzie E, Shivakumar V, Tilton P, Winterbourne SM, Xue S, Kavoussanakis K, Lareau LF, Shah P, Wallace EWJ. 2021. riboviz 2: A flexible and robust ribosome profiling data analysis and visualization workflow. bioRxiv. doi: [10.1101/2021.05.14.443910](https://doi.org/10.1101/2021.05.14.443910).

Wallace, Edward; Anderson, Felicity; Kavoussanakis, Kostas; Jackson, Michael; Shah, Premal; Lareau, Liana; et al. (2021): riboviz: software for analysis and visualization of ribosome profiling datasets. figshare. Software. doi: [10.6084/m9.figshare.12624200](https://doi.org/10.6084/m9.figshare.12624200)

riboviz: analysis and visualization of ribosome profiling datasets, Carja et al., BMC Bioinformatics 2017. doi:[10.1186/s12859-017-1873-8](https://doi.org/10.1186/s12859-017-1873-8).

## Acknowledgements

For contributors and funders, see [Acknowledgements](./docs/acks.md).

For citations of third-party software used by **riboviz**, see [References](./docs/reference/references.md).

## Copyright and License

**riboviz** is Copyright (2016-2021) The University of Edinburgh; Rutgers University; University of California, Berkeley; The University of Pennsylvania.

**riboviz** is released under the [Apache License 2.0](./LICENSE).
