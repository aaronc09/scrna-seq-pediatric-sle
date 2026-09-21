# Review supporting a targeted primary annotation revision

The available evidence supports retaining the defined 93,914-cell T-cell subset in the primary analysis. This is a conservative correction to four previously excluded clusters, not a comprehensive reannotation or experimental validation of all cells. The exact script-14 selection rule is retained without tuning it to statistical significance. Its initial communication results were already known when this further review was undertaken, so adoption as primary is explicitly retrospective.

## Evidence reviewed

Script 16 reads expression for T, helper/naive-memory, CD8, cytotoxic, B, myeloid, platelet, and erythroid marker panels. It compares donor-level and cohort-level QC, stored Scrublet scores, reference predictions, and expression patterns. It does not read communication outcomes. The full tables here include comparison populations and remaining excluded cells, rather than only favorable markers.

| Reviewed population | Cells | Median genes | Median mitochondrial % | Median stored doublet score |
| --- | ---: | ---: | ---: | ---: |
| Pediatric restored CD4 | 44,396 | 965 | 2.10 | 0.064 |
| Pediatric original CD4 | 28,319 | 1,124 | 2.27 | 0.063 |
| Pediatric restored CD8 | 32,084 | 1,077 | 2.64 | 0.059 |
| Pediatric original CD8 | 25,501 | 978 | 2.94 | 0.034 |
| Adult restored CD4 | 13,000 | 982 | 1.96 | 0.068 |
| Adult original CD4 | 10,455 | 1,092 | 2.40 | 0.066 |
| Adult restored CD8 | 4,434 | 1,080 | 2.46 | 0.052 |
| Adult original CD8 | 4,387 | 1,030 | 2.90 | 0.050 |

Both restored CD4 and CD8 populations occur in every one of the 44 pediatric and 12 adult donors. Median mitochondrial percentages within donor/restored-lineage groups range from 0.89% to 4.32%. Thus, these are not cells contributed by only one donor. Detailed donor counts and QC values are in `donor_review.csv`.

All restored cells express CD3D and CD3E by the existing selection rule; this is a criterion, not independent validation. Additional support includes CD8B detection in 87.2% of restored pediatric CD8 and 88.5% of restored adult CD8 cells, compared with 7.4% and 3.6% in restored CD4 cells. IL7R is detected in 74.7% of pediatric and 70.6% of adult restored CD4 cells. CD4 transcript detection is sparse (15.8% and 17.6%), so assignment is supported by reference classes and the broader marker pattern rather than CD4 detection alone. Marker evidence and CellTypist predictions derive from the same RNA data and are not independent assays.

Stored doublet scores overlap those of the retained T-cell populations. Restored pediatric CD8 cells have a higher median than the original CD8 population, which warrants caution; this review does not prove absence of residual doublets. Every reviewed cell already passed the original donor-level Scrublet filter. No new universal doublet-score cutoff was invented, and Scrublet was not rerun on a selected T-cell subset.

## Final rule and its limits

Retain all originally included cells with their original labels. Within pediatric clusters 0/1 and adult clusters 0/4, add cells with CD3D and CD3E detection, an explicitly mapped helper/CD4 or cytotoxic/CD8 CellTypist class, and maximum reference probability >=0.5. Exact class lists and model hash are in the original sensitivity protocol and classifier provenance. Keep RBC/platelet exclusions. Label the other cells in these four clusters as unresolved, not established low-quality cells. No additional ribosomal-percent exclusion is applied.

This targeted rescue applies an additional evidence requirement to previously unannotated cells; it is not a universal QC rule. Only 40–59% of originally retained T cells would meet that same gate. Therefore, using it as a universal quality criterion would be inappropriate. Retained CD8 populations also include reference classes outside the limited rescue mapping (for example, MAIT-like classes); original labels are not claimed to have been comprehensively corrected.

The simple non-T multimarker flags are diagnostic only. Detection of two genes from a lineage panel can reflect shared expression or ambient RNA; these flags occur in roughly 35–42% of restored populations and 40–47% of original T-cell populations. They are not doublet calls and were not used as automatic exclusion rules.

The original four clusters contain 45,420 cells still unresolved after rescue (36,505 pediatric and 8,915 adult). Their exclusion can still introduce selection bias, including bias toward cells with detected CD3 transcripts and higher reference confidence. The primary manuscript must retain this limitation and the original-exclusion sensitivity. No claim is made that the remaining cells are unusable, that the 0.5 cutoff is biologically optimal, or that other annotation schemes would give identical results.

The revised primary results and exact included/excluded labels are under `Results/revised_primary/`. The original-exclusion analysis remains available as a sensitivity comparison.

## Method references

The [CellTypist documentation](https://celltypist.readthedocs.io/en/latest/celltypist.classifier.AnnotationResult.html) describes probability values as sigmoid-transformed decision scores. Here they are used as model confidence, not a calibrated probability that a cell is biologically correct or high quality. The [single-cell best-practices QC chapter](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html) provides the context for reviewing multiple QC measures rather than treating one expression fraction as sufficient evidence of poor quality.
