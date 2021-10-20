# Style guide

This page summarises style guidelines for the riboviz source code, documentation, parameters and files.

* ['riboviz'](#riboviz)
* [Python style](#python-style)
* [R style](#r-style)
* [Nextflow style](#nextflow-style)
* [Configuration parameters](#configuration-parameters)
* [Command-line tools](#command-line-tools)
* [Input and output file names](#input-and-output-file-names)

---

## 'riboviz'

The sofware is called 'riboviz', one word, lower-case, no hyphens.

Please do **not** refer to the software as RiboViz, Riboviz, Ribo-Viz, Ribo-viz or ribo-viz.

---

## Python style

See [Python style](#python-style) in [Developing Python components](./dev-python.md).

---

## R style

See [R style](#r-style) in [Developing R components](./dev-r.md).

---

## Nextflow style

See [Nextflow style](#nextflow-style) in [Developing Nextflow workflow](./dev-nextflow.md).

---

## Configuration parameters

Configuration parameters must be in snake-case i.e., lower-case and delimited by underscores, not hyphens. For example, `adapters` , `asite_disp_length_file`, `orf_gff_file`.

---

## Command-line tools

For consistency with bash and other command-line tools, command-line parameters should be implemented as one, or both, of:

* A single alphanumeric character prefixed by a hyphen. For example, `-v`, `-c 123`, `-s GATC`.
* A kebab-case token i.e., lower-case and delimited by hyphens, not underscores, and prefixed by two hyphens. For example, `--verbose`, `--control=123`, `--match-sequence=GATC`.

---

## Input and output file names

riboviz-specific input and output file names must be in snake-case i.e., lower-case and delimited by underscores, not hyphens. Upper-case is permitted for acronyms e.g., `ORF`, `CDS`, `APE`, `TPMs`, `RNA`.

For example, `vignette_config.yaml`, `read_counts_per_file.tsv`, `TPMs_all_CDS_all_samples.tsv`.
