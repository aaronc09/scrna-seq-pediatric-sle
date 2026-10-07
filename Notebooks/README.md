<!-- reviewer-revision-2026-10-06 -->
The current workflow now ends with **18_REVIEWER_REVISION.py** (isolated LIANA rank correction and fresh inference) and **19_SUMMARIZE_REVIEWER_REVISION.py** (corrected models and reviewer sensitivities). See [the revised output guide](../Results/reviewer_revision/README.txt) for prerequisites and commands. The script-17 publication workflow described below is historical and must be followed by stages 18–19 for corrected rank-based outputs.

# One project code set

This is the sole active project-code directory. Notebook files (`.ipynb`) and Python scripts (`.py`) implement different stages of the same research workflow. They are not duplicate implementations to choose between.

| Stage | Files | Purpose |
| --- | --- | --- |
| Preprocess and establish original annotations | 01A–01F notebooks | QC, original cluster labels, and donor splitting |
| Establish original inference and comparison results | 02A–06B notebooks | Donor LIANA inference, metadata, original primary models, and historical secondary comparisons |
| Original presentation and supporting exports | 07–10 and 12 scripts | Original figures, tables, clinical summaries, and annotation records |
| Recover and verify source annotations | 13 script | Reconcile archived pre-exclusion expression objects with original records |
| Annotation and resource sensitivities | 14–15 scripts | Generate and compare the four donor-level inference scenarios |
| Review and build current primary results | 16–17 scripts | Review restored T cells; fit revised primary and female-only models; generate current figures and tables |
| Shared functions and settings | Unnumbered `.py` modules | Imports used by notebooks and scripts; retain these files |

Follow the root README for prerequisites and run order. Historical original-exclusion stages are retained because they establish comparison results and support reproducibility. They are not a second current primary analysis.

The historical stage-17 workflow is:

```powershell
python Notebooks/16_REVIEW_RESTORED_T_CELLS.py --archive-root "C:/path/to/GSE135779_transfer"
python Notebooks/17_BUILD_REVISED_PRIMARY.py --archive-root "C:/path/to/GSE135779_transfer" --publish
```

Stage 17 with `--publish` produces the earlier annotation-revised outputs; it must be followed by stages 18-19 for the corrected primary analysis. See `Results/revised_primary/README.md` for output locations and prerequisites.

There is deliberately no active script 11: its obsolete submission exporter contained original-analysis values and was removed. Use stages 18-19 for corrected figure/table exports. Local package dependencies under `_local_submission/sensitivity_cache/python_packages/` are third-party libraries, not another set of project code.

Current Table 2 PNG can be regenerated with `python Notebooks/render_table2_png.py`. Numbered supplemental files are packaged by `Notebooks/package_supplemental_tables.py`.
