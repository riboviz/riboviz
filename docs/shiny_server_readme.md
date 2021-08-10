---
title: "Using the shiny visualization tool"
output: 
  html_document:
    df_print: paged
author: "John Favate"
date: "`r Sys.time()`"
---

Riboviz 2.0 has an integrated, but optional data visualization tool written in `R` and using `shiny`. TO use the visualization tool, you can enter the following command at the command line. If you're outside of `riboviz/rscripts`, then you'll need to provide the full path to the `run_shiny_server.R` script.
```{bash, echo = FALSE}
Rscript --vanilla run_shiny_server.R path/to/yaml/configuration/file.yaml
```

Once entered, the command will provide an IP address which you can simply paste into your browsers URL bar and navigate to. If you're running on a headless machine like a server, you will need to first forward the port the shiny server has been created on an be on the same network as the headless machine in order to access it. By default, the port used is `4254` but if you need to give a specific port, you can supply one when entering the command as shown below. This will use the port `1234`.
```{bash, echo = FALSE}
Rscript --vanilla run_shiny_server.R path/to/yaml/configuration/file.yaml 1234
```

# Description of tabs

Globally, the user can choose which samples they would like to view by checking the appropriate boxes at the top of the screen. 

### Read count summary

This tab displays a bar graph showing the number reads over the course of the processing steps.

### TPM summary

This tab display a heatmap with sample to sample pairwise correlations based on $log_{10}(TPM)$ as well as TPM density plots showing the distribution of TPMs for each gene. This tab allows you to enter the name of gene from your data, doing so will draw a vertical black line on the density plots to show you where that gene lies among the distribution of TPMs. 

### Read length distributions

This tab displays read length distributions.

### 3nt periodicity

This tab displays graphs related to the positioning and the frame of the reads.

### Normalized reads per position

This tab displays the normalized read coverage around start and stop sites in each of the samples. 

### Gene specific coverage

This displays read distributions for each codon along a single gene. Reads from each frame are summed to give a single number for each codon and plotted. The user can select the gene and region within that gene to view.

### Ribogrid

This tab displays a meta-feature plot called a ribogrid, The user can select the read lengths and positions they wish to view. Bar graph versions of the ribogrids are shown below. 

### Features

This tab displays scatterplots of TPMs versus other features of the genome such as GC content, transcript length, or other features. The user can select which features to graph against. This plot is only possible when the user has supplied a features file. 
