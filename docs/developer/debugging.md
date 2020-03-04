# Debugging

## Debugging R scripts with appropriate command-line arguments

To debug R scripts such as `generate_stats_figs.R` and `bam_to_h5.R`, they need to be run with the correct command-line arguments to discover the bug. R has good tools for interactive debugging, [explained in Hadley Wickham's chapter on debugging in R](https://adv-r.hadley.nz/debugging.html). However, interactive debugging tools such as `browser()` don't interrupt a call to `Rscript`. Instead you need to modify the call from

```console
Rscript code_to_debug.R --myarg1 value1
```

to

```console
R --args --myarg1 value1
```

then, from the R prompt run

```R
> source('code_to_debug.R')
```

this will accept `debug()` and `browser()` statements run from the interactive R prompt.

For example, in the vignette we call:

```console
Rscript --vanilla rscripts/generate_stats_figs.R --Ncores=1 \
 --MinReadLen=10 --MaxReadLen=50 --Buffer=250 --PrimaryID=Name \
 --dataset=vignette --hdFile=vignette/output/WTnone.h5 \
 --out_prefix=vignette/output/WTnone \
 --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True \
 --dir_out=vignette/output --do_pos_sp_nt_freq=True \
 --t_rna=data/yeast_tRNAs.tsv \
 --codon_pos=data/yeast_codon_pos_i200.RData \
 --features_file=data/yeast_features.tsv \
 --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
 --count_threshold=64 \
 --asite_disp_length_file=vignette/input/asite_disp_length_yeast_standard.txt
```

But to debug a new feature, instead run:

```console
R --vanilla --args --Ncores=1 --MinReadLen=10 --MaxReadLen=50 \
 --Buffer=250 --PrimaryID=Name --dataset=vignette \
 --hdFile=vignette/output/WTnone.h5 --out_prefix=vignette/output/WTnone \
 --orf_fasta=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True \
 --dir_out=vignette/output --do_pos_sp_nt_freq=True \
 --t_rna=data/yeast_tRNAs.tsv \
 --codon_pos=data/yeast_codon_pos_i200.RData \
 --features_file=data/yeast_features.tsv \
 --orf_gff_file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
 --count_threshold=64 \
 --asite_disp_length_file=vignette/input/asite_disp_length_yeast_standard.txt 
```

then 

```R
> source('rscripts/generate_stats_figs.R')
```

For example, to debug a specific line of code, you could add a `browser()` statement in the source first. Alternatively, you could copy and paste the parts of the code you wanted to run, as long as earlier dependencies are run first (packages, importing command arguments, function definitions).
