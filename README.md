# RiboViz

**Ribosome profiling** provides a detailed global snapshot of protein synthesis in a cell.  At its core, this technique makes use of the observation that a translating ribosome protects around 30 nucleotides of the mRNA from nuclease activity.  High-throughput sequencing of these ribosome protected fragments (called ribosome footprints) offers a precise record of the number and location of the ribosomes at the time at which translation is stopped. Mapping the position of the ribosome protected fragments indicates the translated regions within the transcriptome.  Moreover, ribosomes spend different periods of time at different positions, leading to variation in the footprint density along the mRNA transcript. This provides an estimate of how much protein is being produced from each mRNA. Importantly, ribosome profiling is as precise and detailed as RNA sequencing. Even in a short time, since its introduction in 2009, ribosome profiling has been playing a key role in driving biological discovery.

We have developed a bioinformatics tool-kit, **riboviz**, for analyzing and visualizing ribosome profiling datasets. **RiboViz** consists of a comprehensive and flexible backend analysis pipeline and a web application for visualization. The current iteration of **RiboViz* is designed for yeast datasets.

Existing yeast datasets consist of a mix of studies, some of which use elongation inhibitors such as cycloheximide (CHX) and others that flash freeze (FF) the samples to prevent initiation and elongation during sample preparation. In general, each experimental step can potentially introduce biases in processed datasets. **RiboViz** can help identify these biases by allowing users to compare and contrast datasets obtained under different experimental conditions.

All the codes for processing the raw reads is available in this repository.

In addition to the metagenomic analyses, an [R](https://www.r-project.org/)\/[Shiny](https://shiny.rstudio.com/) integration allows the user to select a gene of interest and compare ribosomal densities of its ORF across three data sets.

## Getting started

* [Introduction](./docs/introduction.md)
* [Install prerequisites](./docs/install.md)
* [Quick install scripts](./docs/quick-install.md) (Ubuntu and CentOS only)
* [Map mRNA and ribosome protected reads to transcriptome and collect data intoan HDF5 file](./docs/run-vignette.md). This page describes how you can run a "vignette" of the back-end analysis pipeline, to demostrate RiboViz's capabilities.
* [Prepare data](./docs/prepare-data.md)
* [Structure of HDF5 data](./docs/hdf5-data.md)
