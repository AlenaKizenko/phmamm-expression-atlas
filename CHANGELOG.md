# Changelog

## v1.0.0 (2026-11-19)

First public version.

- Gene search by name, description, eggNOG name, Pfam domain, GO term, KEGG or any identifier; filter by maternal / zygotic class.
- Gene page: expression plot (TPM or counts; mean + replicates, replicate series, or bars; optional log2), similar
  expression profiles (Pearson r), genomic / mRNA / CDS / protein sequences with FASTA download, RefSeq and eggNOG annotation.
- Compare: heatmap of many genes (IDs and/or a search term), coloured by log2(TPM + 1), z-score or TPM, ordered by
  pattern clustering, expression-level clustering, or as entered.
- Maternal / zygotic class (present in the unfertilized egg at >= 1 TPM or not) and peak stage for every nuclear gene.
- Mitochondrial genes (MITOS2) on their own TPM scale.
- Stages table (Hotta staging) and downloads of the full TPM, counts and annotation tables.
- Database: 37,157 genes (37,122 RefSeq + 35 mitochondrial), 15 stages x 2 replicates (ENA PRJEB25149).
