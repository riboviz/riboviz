"""
General FASTA and GFF related constants and functions.
"""

CDS_FEATURE_FORMAT = "{}_CDS"
"""
CDS feature name format for CDS features which do not define ``ID``
or ``Name`` attributes.
"""
START_CODON = "ATG"
""" Canonical start codon. """
STOP_CODONS = ["TAA", "TAG", "TGA"]
""" Canonical stop codons. """
