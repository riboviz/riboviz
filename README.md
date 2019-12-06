# RiboViz

**Ribosome profiling** provides a detailed global snapshot of protein synthesis in a cell.  At its core, this technique makes use of the observation that a translating ribosome protects around 30 nucleotides of the mRNA from nuclease activity.  High-throughput sequencing of these ribosome protected fragments (called ribosome footprints) offers a precise record of the number and location of the ribosomes at the time at which translation is stopped. Mapping the position of the ribosome protected fragments indicates the translated regions within the transcriptome.  Moreover, ribosomes spend different periods of time at different positions, leading to variation in the footprint density along the mRNA transcript. This provides an estimate of how much protein is being produced from each mRNA. Importantly, ribosome profiling is as precise and detailed as RNA sequencing. Even in a short time, since its introduction in 2009, ribosome profiling has been playing a key role in driving biological discovery.

We have developed this bioinformatics toolkit, **RiboViz**, for analyzing and visualizing ribosome profiling datasets. **RiboViz** consists of a comprehensive and flexible analysis pipelin. The current version of **RiboViz** is designed for yeast datasets.

Existing yeast datasets consist of a mix of studies, some of which use elongation inhibitors such as cycloheximide (CHX) and others that flash freeze (FF) the samples to prevent initiation and elongation during sample preparation. In general, each experimental step can potentially introduce biases in processed datasets. **RiboViz** can help identify these biases by allowing users to compare and contrast datasets obtained under different experimental conditions.

All the code for processing the raw reads is available in this repository.

In addition to the metagenomic analyses, an [R](https://www.r-project.org/)\/[Shiny](https://shiny.rstudio.com/) integration provides a web application for visualization which allows the user to select a gene of interest and compare ribosomal densities of its ORF across three data sets.

The current version of **RiboViz** is designed for yeast datasets but it can be customised and used to analyse datasets relating to other organisms.

## Use RiboViz

Introduction:

* [Introduction](./docs/introduction.md)

Quick start:

* [Install prerequisites](./docs/install.md)
* [Quick install scripts](./docs/quick-install.md) (Ubuntu and CentOS only)
* [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./docs/run-vignette.md). Run a "vignette" of the RiboViz workflow to see RiboViz's capabilities.
* [Run UMI extraction, deduplication and demultiplexing examples](./docs/run-dedup-demultiplex-examples.md). Run RiboViz on simulated data, to see how RiboViz handles duplicated and multiplexed data.

Usage:

* [What the RiboViz workflow does](./docs/prep-riboviz-operation.md)
* [Configuring the RiboViz workflow](./docs/prep-riboviz-config.md)
* [Running the RiboViz workflow](./docs/prep-riboviz-running.md)
* [Prepare data](./docs/prepare-data.md)

Reference:

* [create_fastq_simdata.py example FASTQ file generator](./docs/create-fastq-examples.md)
* [demultiplex-fastq.py fastq demultiplexer](./docs/demultiplex-fastq.md)
* [Content and provenance of repository data files](./docs/data.md)
* [Structure of HDF5 data](./docs/hdf5-data.md)

## Develop RiboViz

* [Developer guide](./docs/developer-guide.md)

## Releases

| Release | Description |
| ------- | ----------- |
| [1.1.0](https://github.com/riboviz/RiboViz/releases/tag/1.1.0) | Most recent version prior to commencement of BBSRC/NSF Riboviz project |
| [1.0.0](https://github.com/riboviz/RiboViz/releases/tag/1.0.0) | Associated with Carja et al. (2017) "riboviz: analysis and visualization of ribosome profiling datasets", BMC Bioinformatics, volume 18, article 461 (2017), 25 October 2017, doi: [10.1186/s12859-017-1873-8](https://doi.org/10.1186/s12859-017-1873-8) |
| [0.9.0](https://github.com/riboviz/RiboViz/releases/tag/0.9.0) | Additional code/data associated with the paper below |
| [0.8.0](https://github.com/riboviz/RiboViz/releases/tag/0.8.0) | Associated with Carja et al. (2017) "riboviz: analysis and visualization of ribosome profiling datasets", bioRXiv, 12 January 2017,doi: [10.1101/100032](https://doi.org/10.1101/100032) |

## License

**RiboViz** is released under the [Apache License 2.0](./LICENSE).
