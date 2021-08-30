# Run RiboViz with SGE Executor

This page shows how to run the RiboViz workflow with the SGE executor.

---

## Configuration file

The configuration file (`SGE.config`) provides the configurations and parameters to run the RiboViz workflow on multiple nodes of a HPC system like Eddie.

There is a preset configuration file, `SGE.config`, in the repositoy folder, which allows the RiboViz run on Eddie.
You can also customize the configuration file to allow it run on a different HPC system, or change the preset parameters.
See [Nextflow on Eddie](https://git.ecdf.ed.ac.uk/igmmbioinformatics/nextflow-eddie) for more details on how to create and modify the file.

## Run RiboViz with SGE on Eddie

### Create a interactive session
The workflow management system of RiboViz, Nextflow, uses Java Virtual Machine. The RiboViz can not run on the login node of Eddie, due to the limit of JVM memory limit. The compute nodes of Eddie can not submit jobs to the Eddie's batch system, so multiple-nodes version RiboViz can not run on these nodes too. Therefore, we will use interactive sessions to run multiple-nodes version RiboViz.

To use the interactive session, you should use `qlogin`.

```console
$ qlogin -l h_vmem=16G -pe interactivemem 1
```
The `-pe interactivemem` specifies the parallel environment to be interactive and `1` represents one core is needed.
`-l h_vmem=16G` represents this task needs 16GB memory.

The interactive session are only used to execute the Nextflow and some short tasks, so a one-core environment is enough. However, the memory should be 8GB, 16GB or bigger, because JVM requires much memory.

### Prepare Environment
You should follow the same steps in [run-on-eddie.md](./run-on-eddie.md) to prepare the environment.
You need to downlad the RiboViz, the dataset and all required software.

### Run RiboViz
Running RiboViz with SGE Executor is similar to the instructions in [run-on-eddie.md](./run-on-eddie.md). We only need to add a `-c SGE.config` in the command line parameters.

For example, when we want to run the vignette dataset on multiple nodes, the command is:
```console
nextflow run prep_riboviz.nf -c SGE.config -params-file vignette/vignette_config.yaml -ansi-log false
```
You can see that the only difference is `-c SGE.config`.

### Changable Parameters
Though the multiple-node version RiboViz is designed to be ready to use with out changing any parameters in the existing yaml files, you can still adjust some parameters to achieve higher performance or lower the possibility of failure.
There are a total of three parameters that are changable.

#### default time
You can add a line in the yaml file to control the default time of each tasks in the workflow, for example:
```yaml
default_time: 5h
```
This time will be used as the default value of estimating the time of each task. The default value is 2h (2 hours).

#### time mult
You can add a line in the yaml file to control the time multiplier of each tasks in the workflow, for example:
```yaml
time_mult: 3
```
To guarantee the tasks will not be killed by short of time, the extimation value generated from the time estimation model will be multiplied by this parameter. The default value will not be multiplied by this parameter. The default value is 3.

#### time override
You can add a line in the yaml file to decide to use the time estimation model or use the default value, for example:
```yaml
time_override: false
```
When set to false, the workflow will use a bulit-in time estimation model to estimate the time of each task and fill in the time field of the submission script.
When set to True, the workflow will apply the default time as the estimated time of all tasks.
The default value is false.

## Troubleshooting
### Task get killed by 140 (wall clock limit)
This is because the time field of the submission script of some tasks are too short. You can increase default_time or time_mult or both.

### Long schedule wait time
You may find out that the task are submitted to the batch system, but it keeps waiting in the schedule queue for a long time (i.e. more than 5 hours). The reason is that Eddie is buzy and it can not reserve enough resources for you.

To relieve this, you can:
* Reduce the `num_processes` parameter. It is investigated that most tools in the RiboViz workflow can achieve a high speedup at 4-8 processes.
* Decrease the default_time or time_mult or both. This will make the time field in the submission script shorter, so the tasks can fill in more empty slots. However, this may cause tasks get killed before it finishes.
