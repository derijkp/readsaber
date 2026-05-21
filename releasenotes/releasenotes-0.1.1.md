
Changes/features
- simplifyschema - keep forward (iso reverse)
- in case of multiple overlapping matches in -refseq, pick the one with the highest mapquality: (useful when using -refseqannot '$chromosome')

Fixes
- fixed error in handling of overlapping annotations
- readsaber_graph - add gray if not enough colors for all elements
- readsaber_graph - get elements from files so extra elements from using -refseqannot '$chromosome' are not missed
