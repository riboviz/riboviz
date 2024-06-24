# Check FASTA and GFF files for coding sequence (CDS) features

`check_fasta_gff` checks FASTA and GGF files for coding sequence (CDS) features. It can be run as follows:

```console
$ check_fasta_gff [-h] \
         -f FASTA -g GFF [-o FEATURES_ISSUES] \
         [--use-feature-name] \
         [--feature-format FEATURE_FORMAT]
         [--start-codon START_CODON [START_CODON ...]] \
	 [-v]
```

where:

* `-f FASTA`, `--fasta FASTA`: FASTA file input.
* `-g GFF`, `--gff GFF`: GFF3 file input.
* `-o FEATURES_ISSUES`, `--features-issues FEATURES_ISSUES`: Issues file output (default `features_issues.tsv`).
* `--use-feature-name`: If a CDS feature defines both `ID` and `Name` attributes then use `Name` in reporting, otherwise use `ID` (default `false`).
* `--feature-format FEATURE_FORMAT`: Feature name format for features which do not define `ID` or `Name` attributes. This format is applied to the sequence ID to create a feature name (default `{}_CDS`).
* `--start-codon START_CODON [START_CODON ...]`: Allowable start codons (default `ATG`).
* `-v`: Print information on each issue (default `false`)

Issues are both reported to the console and saved in an issues file.

The issues file is a file of tab-separated values. The file has columns:

* `Sequence`: sequence ID. This value is one of:
  - Value of sequence ID.
  - `*`, if the issue relates to multiple sequences.
* `Feature`: feature ID.
  - Empty, if the issue relates to the sequence, not the feature.
  - `*`, if the issue relates to multiple features.          
  - Value of `ID` attribute for feature, if defined.
  - Value of `Name` attribute for feature, if defined, and if `ID` is undefined.
  - Sequence ID formatted using `feature_format` (default `{}_CDS`) if both `ID` and `Name` are undefined.
* `Issue`: issue type. See below.
* `Data`: issue data or empty, see below.

The following issue types are reported for every CDS annotated in the GFF:

* `IncompleteFeature`: The CDS has a length not divisible by 3.
* `NoStartCodon` : The CDS does not start with a start codon (`ATG` or, those provided via `--start_codon`). The supplementary issue data is the actual codon found.
* `NoStopCodon` : The CDS does not end with a stop codon  (`TAG`, `TGA`, `TAA`). The supplementary issue data is the actual codon found.
* `InternalStopCodon`: The CDS has internal stop codons.
* `NoIdName`: The CDS has no `ID` or `Name` attribute.
* `DuplicateFeatureId`: The CDS has a non-unique `ID` attribute (attributes are expected to be unique within the scope of a GFF file).
* `DuplicateFeatureIds`: Related to the above, multiple CDSs have non-unique `ID` attributes. This summarises the count of all CDSs that share a common `ID` attribute. For this issue, the sequence `ID` attribute is `*`. The supplementary issue data is a count of the number of features with the same ID.

The following issues are reported for sequences defined in the GFF file:

* `MultipleCDS`: The sequence has multiple CDS. For this issue, the feature `ID` attribute is `*`. The supplementary issue data is a count of the number of CDSs found.
* `SequenceNotInFASTA` : The sequence has a feature in the GFF file but the sequence is not in the FASTA file. For this issue, the feature `ID` attribute is rmpty.
* `SequenceNotInGFF`: The sequence is in the FASTA file but has no features in the GFF file. For this issue, the feature `ID` attribute is empty.

Example:

```console
$ check_fasta_gff \
    -f vignette/input/yeast_YAL_CDS_w_250utrs.fa \
    -g vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
    -o check_vignette_YAL.tsv
...
Configuration:
fasta_file	vignette/input/yeast_YAL_CDS_w_250utrs.fa
gff_file	vignette/input/yeast_YAL_CDS_w_250utrs.gff3
start_codons	['ATG']

Metadata:
NumSequences:	68
NumFeatures:	204
NumCDSFeatures:	68

Issue summary:
Issue	Count
NoStopCodon	1
InternalStopCodon	1
NoStartCodon	1
DuplicateFeatureId	0
SequenceNotInFASTA	0
MultipleCDS	0
DuplicateFeatureIds	0
IncompleteFeature	0
NoIdName	0
SequenceNotInGFF	0

$ cat check_vignette_YAL.tsv
...
# fasta_file: vignette/input/yeast_YAL_CDS_w_250utrs.fa
# gff_file: vignette/input/yeast_YAL_CDS_w_250utrs.gff3
# start_codons: ['ATG']
# NumSequences: 68
# NumFeatures: 204
# NumCDSFeatures: 68
# NoStopCodon: 1
# InternalStopCodon: 1
# NoStartCodon: 1
# DuplicateFeatureId: 0
# SequenceNotInFASTA: 0
# MultipleCDS: 0
# DuplicateFeatureIds: 0
# IncompleteFeature: 0
# NoIdName: 0
# SequenceNotInGFF: 0
Sequence	Feature	Issue	Data
YAL001C	YAL001C	NoStartCodon	AAA
YAL001C	YAL001C	NoStopCodon	TTT
YAL001C	YAL001C	InternalStopCodon	
```

Example with `-v` verbose mode:

```console
$ check_fasta_gff \
    -f vignette/input/yeast_YAL_CDS_w_250utrs.fa \
    -g vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
    -o check_vignette_YAL.tsv -v
...
Configuration:
fasta_file	vignette/input/yeast_YAL_CDS_w_250utrs.fa
gff_file	vignette/input/yeast_YAL_CDS_w_250utrs.gff3
start_codons	['ATG']

Metadata:
NumSequences:	68
NumFeatures:	204
NumCDSFeatures:	68

Issue summary:
Issue	Count
InternalStopCodon	1
NoStopCodon	1
NoStartCodon	1
DuplicateFeatureId	0
SequenceNotInGFF	0
IncompleteFeature	0
MultipleCDS	0
SequenceNotInFASTA	0
DuplicateFeatureIds	0
NoIdName	0

Issue details:
Sequence YAL001C feature YAL001C doesn't start with a recognised start codon but with AAA
Sequence YAL001C feature YAL001C doesn't end with a recognised stop codon but with TTT
Sequence YAL001C feature YAL001C has an internal stop codon
```

Example with `--start-codon`:

```console
$ check_fasta_gff \
     -f vignette/input/yeast_YAL_CDS_w_250utrs.fa \
     -g vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
     -o check_vignette_YAL.tsv \
     --start-codon ATG AAA -v
...
Configuration:
fasta_file	vignette/input/yeast_YAL_CDS_w_250utrs.fa
gff_file	vignette/input/yeast_YAL_CDS_w_250utrs.gff3
start_codons	['ATG', 'AAA']

Metadata:
NumSequences:	68
NumFeatures:	204
NumCDSFeatures:	68

Issue summary:
Issue	Count
NoStopCodon	1
InternalStopCodon	1
DuplicateFeatureId	0
NoStartCodon	0
SequenceNotInFASTA	0
IncompleteFeature	0
MultipleCDS	0
DuplicateFeatureIds	0
NoIdName	0
SequenceNotInGFF	0

Issue details:
Sequence YAL001C feature YAL001C doesn't end with a recognised stop codon but with TTT
Sequence YAL001C feature YAL001C has an internal stop codon
```

An example with more issues:

```console
$ check_fasta_gff \
     -f data/yeast_CDS_w_250utrs.fa \
     -g data/yeast_CDS_w_250utrs.gff3 \
     -o check_data_CDS.tsv -v
...
Configuration:
fasta_file	data/yeast_CDS_w_250utrs.fa
gff_file	data/yeast_CDS_w_250utrs.gff3
start_codons	['ATG']

Metadata:
NumSequences:	5812
NumFeatures:	17436
NumCDSFeatures:	5812

Issue summary:
Issue	Count
InternalStopCodon	17
NoStartCodon	1
NoIdName	0
SequenceNotInGFF	0
MultipleCDS	0
SequenceNotInFASTA	0
DuplicateFeatureId	0
IncompleteFeature	0
NoStopCodon	0
DuplicateFeatureIds	0

Issue details:
Sequence Q0050 feature Q0050 has an internal stop codon
Sequence Q0055 feature Q0055 has an internal stop codon
Sequence Q0060 feature Q0060 has an internal stop codon
Sequence Q0065 feature Q0065 has an internal stop codon
Sequence Q0070 feature Q0070 has an internal stop codon
Sequence Q0045 feature Q0045 has an internal stop codon
Sequence Q0075 feature Q0075 doesn't start with a recognised start codon but with ATA
Sequence Q0075 feature Q0075 has an internal stop codon
Sequence Q0085 feature Q0085 has an internal stop codon
Sequence Q0110 feature Q0110 has an internal stop codon
Sequence Q0115 feature Q0115 has an internal stop codon
Sequence Q0120 feature Q0120 has an internal stop codon
Sequence Q0105 feature Q0105 has an internal stop codon
Sequence Q0140 feature Q0140 has an internal stop codon
Sequence Q0160 feature Q0160 has an internal stop codon
Sequence Q0250 feature Q0250 has an internal stop codon
Sequence Q0255 feature Q0255 has an internal stop codon
Sequence Q0275 feature Q0275 has an internal stop codon
```
```console
$ cat check_data_CDS.tsv
...
# fasta_file: data/yeast_CDS_w_250utrs.fa
# gff_file: data/yeast_CDS_w_250utrs.gff3
# start_codons: ['ATG']
# NumSequences: 5812
# NumFeatures: 17436
# NumCDSFeatures: 5812
# InternalStopCodon: 17
# NoStartCodon: 1
# MultipleCDS: 0
# SequenceNotInGFF: 0
# IncompleteFeature: 0
# NoStopCodon: 0
# DuplicateFeatureId: 0
# NoIdName: 0
# DuplicateFeatureIds: 0
# SequenceNotInFASTA: 0
Sequence	Feature	Issue	Data
Q0050	Q0050	InternalStopCodon	
Q0055	Q0055	InternalStopCodon	
Q0060	Q0060	InternalStopCodon	
Q0065	Q0065	InternalStopCodon	
Q0070	Q0070	InternalStopCodon	
Q0045	Q0045	InternalStopCodon	
Q0075	Q0075	NoStartCodon	ATA
Q0075	Q0075	InternalStopCodon	
Q0085	Q0085	InternalStopCodon	
Q0110	Q0110	InternalStopCodon	
Q0115	Q0115	InternalStopCodon	
Q0120	Q0120	InternalStopCodon	
Q0105	Q0105	InternalStopCodon	
Q0140	Q0140	InternalStopCodon	
Q0160	Q0160	InternalStopCodon	
Q0250	Q0250	InternalStopCodon	
Q0255	Q0255	InternalStopCodon	
Q0275	Q0275	InternalStopCodon
```
