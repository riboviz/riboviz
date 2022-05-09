# Create job submission script from template

`create_job_script` has a script that can take a template job submission script and customise it with values provided via the command-line or within a riboviz YAML configuration file. The script replaces configuration files, delimited by `%%`-delimited tokens (for example `%%job_num_cpus%%`) with configuration values.

The following configuration parameters can be specified within the template:

| Parameter | Description | Mandatory | Default |
| --------- | ----------- | --------- | ------- |
| `job-email-events` | Events triggering emails about batch job. Any combination of `b`(begin), `e` (end), `a` (abort), `s` (suspend). | No | `beas` |
| `job-email` | E-mail address for batch job events | No | `null` |
| `job-memory` | Requested memory for batch job | No | `8G` |
| `job-name` | Name of batch job | No | `riboviz` |
| `job-num-cpus` | Requested number of CPUs for batch job | No | `4` |
| `job-parallel-env` | Requested parallel environment for batch job. One of `sharedmem`, `mpi`, `scatter`, `gpu`, See University of Edinburgh ECDF Linux Compute Cluster, [parallel environments](https://www.wiki.ed.ac.uk/display/ResearchServices/Parallel+Environments) | No | `mpi` |
| `job-runtime` | Maximum runtime for batch job | No | `48:00:00` |
| `nextflow-dag-file` | Nextflow DAG file | No | `nextflow-dag.html` |
| `nextflow-report-file` | Nextflow report file | No | `nextflow-report.html` |
| `nextflow-timeline-file` | Nextflow timeline file | No | `nextflow-timeline.html` |
| `nextflow-trace-file` | Nextflow trace file | No | `nextflow-trace.tsv` |
| `nextflow-work-dir` | Nextflow work directory | No | `work` |
| `validate-only ` | Validate configuration, check that mandatory parameters have been provided and that input files exist, then exit without running the workflow? | No | `false` |

In addition the following parameters can be specified within the template:

| Parameter | Description |
| --------- | ----------- |
| `config-file` | riboviz YAML configuration file |
| `r-libs` | Location of R libraries required by riboviz |

As an example, a template job submission script for [Eddie](https://www.ed.ac.uk/information-services/research-support/research-computing/ecdf/high-performance-computing), The University of Edinburgh ECDF Linux Compute Cluster, is in `jobs/eddie-template.sh`.

`create_job_script` can be used to customise a template as follows:

```console
$ create_job_script [-h] \
  -i [INPUT_FILE] [-o [OUTPUT_FILE]] [-t [TOKEN_TAG]] \
  --r-libs R_LIBS --config-file [CONFIG_FILE] \
  [--job-email JOB_EMAIL] \
  [--job-email-events JOB_EMAIL_EVENTS] \
  [--job-memory JOB_MEMORY] \
  [--job-name JOB_NAME] \
  [--job-num-cpus JOB_NUM_CPUS] \
  [--job-runtime JOB_RUNTIME] \
  [--nextflow-dag-file NEXTFLOW_DAG_FILE] \
  [--nextflow-report-file NEXTFLOW_REPORT_FILE] \
  [--nextflow-resume] \
  [--nextflow-timeline-file NEXTFLOW_TIMELINE_FILE] \
  [--nextflow-trace-file NEXTFLOW_TRACE_FILE]
  [--nextflow-work-dir NEXTFLOW_WORK_DIR] \
  [--validate-only]
```

where:

* `-i INPUT_FILE`, `--input_file INPUT_FILE`: Job submission script template.
* `-o [OUTPUT_FILE]`, `--output [OUTPUT_FILE]`: Job submission script (if not provided then the job submission script is printed to standard output).
* `-t [TOKEN_TAG]`, `--token-tag [TOKEN_TAG]`: Tag marking up tokens for replacement in job submission script template. Default `%%`.
* `--config-file [CONFIG_FILE]`: YAML configuration file.
* `r_libs R_LIBS`: Location of R libraries required by riboviz.

Default values, defined in code, are overriden by any values provided in the YAML configuration file. These, in turn, are overridden by any values provided via the command-line.

Examples:

```console
$ create_job_script -i jobs/eddie-template.sh \
    -o run_B-Sc_2012.sh \
    --config-file ~/data-folder/Brar_2012_Meiosis_RPF_6-samples_CDS_w_250utrs_config.yaml \
    --r-libs /exports/csce/eddie/biology/groups/wallace_rna/Rlibrary \
    --job-name B-Sc_2012 --job-runtime 48:00:00 \
    --job-memory 8G --job-num-cpus 16 \
    --validate-only
```
```console
$ create_job_script -i jobs/eddie-template.sh \
    -o run_B-Sc_2012.sh \
    --config-file ~/data-folder/Brar_2012_Meiosis_RPF_6-samples_CDS_w_250utrs_config.yaml \
    --r-libs /exports/csce/eddie/biology/groups/wallace_rna/Rlibrary \
    --job-name B-Sc_2012 --job-runtime 48:00:00 \
    --job-memory 8G --job-num-cpus 16
```
```console
$ create_job_script -i jobs/eddie-template.sh \
    -o run_W-Cn-H99_2020.sh \
    --config-file ~/data-folder/Brar_2012_Meiosis_RPF_6-samples_CDS_w_250utrs_config.yaml \
    --r-libs /exports/csce/eddie/biology/groups/wallace_rna/Rlibrary \
    --job-name W-Cn-H99_2020 --job-runtime 48:00:00 \
    --job-memory 8G --job-num-cpus 16 \
    --validate-only
```
```console
$ create_job_script -i jobs/eddie-template.sh \
    -o run_W-Cn-H99_2020.sh \
    --config-file ~/data-folder/Brar_2012_Meiosis_RPF_6-samples_CDS_w_250utrs_config.yaml \
    --r-libs /exports/csce/eddie/biology/groups/wallace_rna/Rlibrary \
    --job-name W-Cn-H99_2020 --job-runtime 48:00:00 \
    --job-memory 8G --job-num-cpus 16
```
