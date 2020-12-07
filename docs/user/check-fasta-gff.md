# Check FASTA and GFF files for compatibility

`riboviz.tools.check_fasta_gff` checks FASTA and GGF files for equality

```console
$ python -m riboviz.tools.check_fasta_gff [-h] \
         -f FASTA -g GFF [-o FEATURES_ISSUES] \
         [--feature-format FEATURE_FORMAT]
```

where:

* `-f FASTA`, `--fasta FASTA`: FASTA file input.
* `-g GFF`, `--gff GFF`: GFF3 file input.
* `-o FEATURES_ISSUES`, `--features-issues FEATURES_ISSUES`: Issues file output (default `features_issues.tsv`).
* `--feature-format FEATURE_FORMAT`: Feature name format for features which do not define `ID` or `Name` attributes. This format is applied to the sequence ID to create a feature name (default `{}_CDS`).

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

The following issue types are reported for CDSs defined in the GFF file:

* `IncompleteFeature`: The CDS has a length not divisible by 3.
* `NoATGStartCodon` : The CDS does not start with a start codon (`ATG`).
* `NoStopCodon` : The CDS does not end with a stop codon (`TAG`, `TGA`, `TAA`).
* `InternalStopCodon`: The CDS has internal stop codons.
* `NoIdName`: The CDS has no `ID` or `Name` attribute.
* `DuplicateFeatureId`: The CDS has a non-unique `ID` attribute (attributes are expected to be unique within the scope of a GFF file).
* `DuplicateFeatureIds`: Related to the above, multiple CDSs have non-unique `ID` attributes. This summarises the count of all CDSs that share a common `ID` attribute. For this issue, the sequence `ID` attribute is `*`. The supplementary issue data data is a count of the number of features with the same ID.

The following issues are reported for sequences defined in the GFF file:

* `MultipleCDS`: The sequence has multiple CDS. For this issue, the feature `ID` attribute is `*`.
* `SequenceNotInFASTA` : The sequence has a feature in the GFF file but the sequence is not in the FASTA file. For this issue, the feature `ID` attribute is rmpty.
* `SequenceNotInGFF`: The sequence is in the FASTA file but has no features in the GFF file. For this issue, the feature `ID` attribute is empty.

Examples:

```console
$ python -m riboviz.tools.check_fasta_gff \
    -f vignette/input/yeast_YAL_CDS_w_250utrs.fa \
    -g vignette/input/yeast_YAL_CDS_w_250utrs.gff3 \
     -o check_vignette_YAL.tsv
...
Sequence YAL001C feature YAL001C doesn't start with ATG
Sequence YAL001C feature YAL001C doesn't stop at end
Sequence YAL001C feature YAL001C has internal STOP
$ cat check_vignette_YAL.tsv 
...
Sequence	Feature	Issue	Data
YAL001C	YAL001C	NoATGStartCodon	
YAL001C	YAL001C	NoStopCodon	
YAL001C	YAL001C	InternalStopCodon	
```
```console
$ python -m riboviz.tools.check_fasta_gff \
    -f data/yeast_CDS_w_250utrs.fa \
    -g data/yeast_CDS_w_250utrs.gff3 \
    -o check_data_CDS.tsv
Sequence Q0050 feature Q0050 has internal STOP
Sequence Q0055 feature Q0055 has internal STOP
Sequence Q0060 feature Q0060 has internal STOP
Sequence Q0065 feature Q0065 has internal STOP
Sequence Q0070 feature Q0070 has internal STOP
Sequence Q0045 feature Q0045 has internal STOP
Sequence Q0075 feature Q0075 doesn't start with ATG
Sequence Q0075 feature Q0075 has internal STOP
Sequence Q0085 feature Q0085 has internal STOP
Sequence Q0110 feature Q0110 has internal STOP
Sequence Q0115 feature Q0115 has internal STOP
Sequence Q0120 feature Q0120 has internal STOP
Sequence Q0105 feature Q0105 has internal STOP
Sequence Q0140 feature Q0140 has internal STOP
Sequence Q0160 feature Q0160 has internal STOP
Sequence Q0250 feature Q0250 has internal STOP
Sequence Q0255 feature Q0255 has internal STOP
Sequence Q0275 feature Q0275 has internal STOP
```
