"""Generate reproducible manuscript and supplementary summary tables."""

from pathlib import Path
import textwrap

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from publication_utils import FIGURE_DPI, save_publication_figure


BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS = BASE_DIR / "Results"
QC = RESULTS / "qc"
SEVERITY = RESULTS / "severity_analysis"
MANUSCRIPT = RESULTS / "manuscript_tables"
MANUSCRIPT.mkdir(parents=True, exist_ok=True)


GROUP_LABELS = {
    "cHD": "Pediatric healthy",
    "cSLE": "Pediatric SLE",
    "aHD": "Adult healthy",
    "aSLE": "Adult SLE",
}


def _wrap_table_value(value: object, width: int = 24) -> str:
    """Wrap long table entries without changing the underlying CSV data."""
    if pd.isna(value):
        return "NA"
    return "\n".join(textwrap.wrap(str(value), width=width)) or ""


def save_table_png(
    table: pd.DataFrame,
    output_stem: Path,
    *,
    title: str,
    rows_per_page: int = 24,
    wrap_width: int = 24,
    large_text: bool = False,
    embed_title: bool = True,
) -> list[Path]:
    """Render a DataFrame as one or more publication-ready 600-DPI PNG files."""
    if table.empty:
        raise ValueError(f"Cannot render an empty table: {title}")

    # Remove only older pages owned by this exact table export.
    for old_page in output_stem.parent.glob(f"{output_stem.name}*.png"):
        old_page.unlink()

    pages = []
    page_count = int(np.ceil(len(table) / rows_per_page))
    for page_index, start in enumerate(range(0, len(table), rows_per_page), start=1):
        page = table.iloc[start : start + rows_per_page].copy()
        display = page.map(lambda value: _wrap_table_value(value, wrap_width))
        max_lines = max(
            1,
            max(str(value).count("\n") + 1 for value in display.to_numpy().ravel()),
        )
        if large_text:
            figure_height = max(
                4.5,
                2.0 + len(display) * (0.55 + 0.22 * (max_lines - 1)),
            )
            figure_width = 12.0
            title_fontsize = 15
            # The original eight-column manuscript tables used 6.5-point text.
            # Short display labels and wrapping permit substantially larger type.
            table_fontsize = 13
        else:
            figure_height = max(
                2.8,
                1.25 + len(display) * (0.32 + 0.13 * (max_lines - 1)),
            )
            figure_width = max(10.0, min(20.0, 1.65 * len(display.columns)))
            title_fontsize = 12
            table_fontsize = 7.5 if len(display.columns) <= 6 else 6.5
        fig, ax = plt.subplots(figsize=(figure_width, figure_height))
        ax.axis("off")
        page_title = title if page_count == 1 else f"{title} (part {page_index} of {page_count})"
        if embed_title:
            ax.set_title(
                page_title,
                fontsize=title_fontsize,
                fontweight="bold",
                pad=12,
                loc="left",
            )
        display_labels = {
            "group": "Group",
            "n_donors": "Donors",
            "age_median_iqr": "Age, median (IQR)",
            "female_male": "Female/male",
            "batches_n": "Batches (n)",
            "retained_cells": "Retained cells",
            "sledai_available_n": "SLEDAI available",
            "sledai_median_iqr": "SLEDAI, median (IQR)",
        }
        artist = ax.table(
            cellText=display.values,
            colLabels=[
                _wrap_table_value(display_labels.get(column, column), 16)
                for column in display.columns
            ],
            cellLoc="left",
            colLoc="left",
            loc="upper left",
            bbox=[0, 0, 1, 0.94 if embed_title else 1],
        )
        artist.auto_set_font_size(False)
        artist.set_fontsize(table_fontsize)
        for (row, _column), cell in artist.get_celld().items():
            cell.set_edgecolor("#B8BEC6")
            cell.set_linewidth(0.45)
            if row == 0:
                cell.set_facecolor("#DCE6F1")
                cell.set_text_props(weight="bold")
            elif row % 2 == 0:
                cell.set_facecolor("#F5F7FA")
        filename = (
            output_stem.with_suffix(".png")
            if page_count == 1
            else output_stem.with_name(f"{output_stem.name}_part{page_index:02d}").with_suffix(".png")
        )
        save_publication_figure(filename, fig)
        plt.close(fig)
        pages.append(filename)
    return pages


def cohort_table() -> None:
    characteristics = pd.read_csv(QC / "cohort_characteristics_summary.csv")
    confounders = pd.read_csv(QC / "potential_confounders_by_group.csv")
    metadata = pd.read_csv(SEVERITY / "patient_clinical_metadata.csv")
    counts = pd.concat(
        [
            pd.read_csv(QC / "child_sample_celltype_counts.csv"),
            pd.read_csv(QC / "adult_sample_celltype_counts.csv"),
        ],
        ignore_index=True,
    )
    retained = (
        metadata[["sample", "Groups"]]
        .merge(counts, on="sample", validate="one_to_many")
        .groupby("Groups")["n_cells"]
        .sum()
    )
    rows = []
    for group in ["cHD", "cSLE", "aHD", "aSLE"]:
        values = characteristics.loc[characteristics["group"].eq(group)].iloc[0]
        batches = confounders.loc[
            confounders["Groups"].eq(group) & confounders["field"].eq("Batch")
        ]
        rows.append(
            {
                "group": group,
                "n_donors": int(values["n_donors"]),
                "age_median_iqr": (
                    f'{values["age_median"]:g} '
                    f'({values["age_q1"]:g}-{values["age_q3"]:g})'
                ),
                "female_male": f'{int(values["female_n"])}/{int(values["male_n"])}',
                "batches_n": "; ".join(
                    f'{row["level"]}:{int(row["n"])}' for _, row in batches.iterrows()
                ),
                "retained_cells": int(retained[group]),
                "sledai_available_n": int(values["sledai_available_n"]),
                "sledai_median_iqr": (
                    f'{values["sledai_median"]:g} '
                    f'({values["sledai_q1"]:g}-{values["sledai_q3"]:g})'
                    if pd.notna(values["sledai_median"])
                    else "NA"
                ),
            }
        )
    table = pd.DataFrame(rows)
    table["group"] = table["group"].map(GROUP_LABELS)
    table.to_csv(MANUSCRIPT / "table1_cohort_characteristics.csv", index=False)


def cell_type_coverage() -> None:
    metadata = pd.read_csv(SEVERITY / "patient_clinical_metadata.csv")
    counts = pd.concat(
        [
            pd.read_csv(QC / "child_sample_celltype_counts.csv"),
            pd.read_csv(QC / "adult_sample_celltype_counts.csv"),
        ],
        ignore_index=True,
    )
    merged = metadata[["sample", "Groups"]].merge(counts, on="sample")
    rows = []
    for (group, cell_type), values in merged.groupby(["Groups", "cell_type"]):
        rows.append(
            {
                "group": group,
                "cell_type": cell_type,
                "total_donors": metadata.loc[metadata["Groups"].eq(group), "sample"].nunique(),
                "donors_with_any_cells": values["sample"].nunique(),
                "donors_with_at_least_10_cells": int((values["n_cells"] >= 10).sum()),
                "minimum_cells_among_donors_with_any": int(values["n_cells"].min()),
            }
        )
    pd.DataFrame(rows).to_csv(MANUSCRIPT / "supplementary_cell_type_coverage.csv", index=False)


def significant_interaction_summaries() -> None:
    significant = pd.read_csv(SEVERITY / "primary_significant_interactions_fdr05.csv")
    significant["source"].value_counts().rename_axis("source").reset_index(
        name="n_significant_interactions"
    ).to_csv(MANUSCRIPT / "significant_sources.csv", index=False)
    significant["target"].value_counts().rename_axis("target").reset_index(
        name="n_significant_interactions"
    ).to_csv(MANUSCRIPT / "significant_targets.csv", index=False)
    significant.groupby(["source", "target"]).size().sort_values(ascending=False).rename(
        "n_significant_interactions"
    ).reset_index().to_csv(MANUSCRIPT / "significant_source_target_pairs.csv", index=False)

    families = significant.assign(
        ligand_receptor=(
            significant["ligand_complex"].astype(str)
            + "-"
            + significant["receptor_complex"].astype(str)
        )
    )
    families.groupby("ligand_receptor").agg(
        n_significant_interactions=("interaction_id", "size"),
        n_positive=("interaction_coef", lambda values: int((values > 0).sum())),
        n_negative=("interaction_coef", lambda values: int((values < 0).sum())),
    ).sort_values("n_significant_interactions", ascending=False).reset_index().to_csv(
        MANUSCRIPT / "significant_ligand_receptor_families.csv", index=False
    )

    score_frames = []
    for group in ["cH", "cSLE", "aH", "aSLE"]:
        scores = pd.read_csv(SEVERITY / f"{group}_per_patient_scores.csv")
        scores["group"] = group
        score_frames.append(scores)
    scores = pd.concat(score_frames, ignore_index=True)
    group_means = (
        scores.loc[scores["interaction_id"].isin(significant["interaction_id"])]
        .groupby(["interaction_id", "group"])["score"]
        .agg(["mean", "count"])
        .reset_index()
    )
    means_wide = group_means.pivot(index="interaction_id", columns="group", values="mean")
    counts_wide = group_means.pivot(index="interaction_id", columns="group", values="count")
    for group in ["cH", "cSLE", "aH", "aSLE"]:
        means_wide[f"n_{group}"] = counts_wide[group]
    means_wide["pediatric_sle_minus_healthy"] = means_wide["cSLE"] - means_wide["cH"]
    means_wide["adult_sle_minus_healthy"] = means_wide["aSLE"] - means_wide["aH"]
    means_wide["difference_in_differences_from_means"] = (
        means_wide["pediatric_sle_minus_healthy"]
        - means_wide["adult_sle_minus_healthy"]
    )
    count_columns = ["n_cH", "n_cSLE", "n_aH", "n_aSLE"]
    supplementary = significant.merge(
        means_wide[count_columns].reset_index(), on="interaction_id", validate="one_to_one"
    )
    if not supplementary[count_columns].sum(axis=1).eq(supplementary["n_observations"]).all():
        raise ValueError("Per-group donor counts do not match model observation counts")
    supplementary.assign(
        direction=np.where(
            supplementary["interaction_coef"] > 0,
            "positive (more positive pediatric SLE-associated change)",
            "negative (more negative pediatric SLE-associated change)",
        )
    ).to_csv(MANUSCRIPT / "supplementary_all_fdr05_interactions.csv", index=False)
    significant[
        [
            "interaction_id",
            "source",
            "target",
            "ligand_complex",
            "receptor_complex",
            "interaction_coef",
            "fdr_q_value",
        ]
    ].merge(means_wide.reset_index(), on="interaction_id", how="left").to_csv(
        MANUSCRIPT / "significant_interaction_group_means.csv", index=False
    )

    representative_keys = [
        ("Classical Monocytes", "NK Cells", "TIMP2", "ITGB1"),
        ("CD4 T Cells", "NK Cells", "CD48", "CD244"),
        ("Non-classical Monocytes", "CD4 T Cells", "SELPLG", "ITGB2"),
        ("Classical Monocytes", "NK Cells", "TGFB1", "ITGB1"),
        ("Non-classical Monocytes", "CD4 T Cells", "CD48", "CD2"),
        ("B Cells", "Non-classical Monocytes", "CD52", "SIGLEC10"),
        ("CD8 T Cells", "Non-classical Monocytes", "CD52", "SIGLEC10"),
    ]
    merged = significant.merge(means_wide.reset_index(), on="interaction_id", how="left")
    representative_rows = []
    for source, target, ligand, receptor in representative_keys:
        match = merged.loc[
            merged["source"].eq(source)
            & merged["target"].eq(target)
            & merged["ligand_complex"].eq(ligand)
            & merged["receptor_complex"].eq(receptor)
        ]
        if len(match) != 1:
            raise ValueError(
                f"Expected one representative result for {source} -> {target}, "
                f"{ligand}-{receptor}; found {len(match)}."
            )
        row = match.iloc[0]
        pediatric_change = row["cSLE"] - row["cH"]
        adult_change = row["aSLE"] - row["aH"]
        if pediatric_change > 0 and adult_change < 0:
            pattern = "Increased in children but decreased in adults"
        elif pediatric_change < 0 and adult_change > 0:
            pattern = "Decreased in children but increased in adults"
        elif pediatric_change < 0 and abs(pediatric_change) > abs(adult_change):
            pattern = "Decreased mainly in children"
        else:
            pattern = "Different SLE-associated changes between cohorts"
        representative_rows.append(
            {
                "Source to target": f"{source} to {target}",
                "Ligand-receptor": f"{ligand}-{receptor}",
                "Coefficient": f'{row["interaction_coef"]:.3f}',
                "95% CI": (
                    f'{row["interaction_ci95_low"]:.3f} to '
                    f'{row["interaction_ci95_high"]:.3f}'
                ),
                "FDR": f'{row["fdr_q_value"]:.5f}',
                "Pediatric healthy to SLE": f'{row["cH"]:.3f} to {row["cSLE"]:.3f}',
                "Adult healthy to SLE": f'{row["aH"]:.3f} to {row["aSLE"]:.3f}',
                "Plain-language pattern": pattern,
            }
        )
    pd.DataFrame(representative_rows).to_csv(
        MANUSCRIPT / "table2_representative_interactions.csv", index=False
    )


def adult_exclusions() -> None:
    pd.DataFrame(
        [
            {
                "original_study_id": "aHD2",
                "original_group": "adult healthy",
                "reason_not_analyzed": (
                    "No corresponding unambiguous GEO GSM/raw matrix mapping was present "
                    "in the transferred project files."
                ),
            },
            {
                "original_study_id": "aSLE8",
                "original_group": "adult SLE",
                "reason_not_analyzed": (
                    "No corresponding unambiguous GEO GSM/raw matrix mapping was present "
                    "in the transferred project files."
                ),
            },
        ]
    ).to_csv(MANUSCRIPT / "supplementary_adult_sample_exclusions.csv", index=False)


def severity_summary() -> None:
    rows = []
    for analysis, filename in [
        ("pediatric SLEDAI correlation", "cSLE_severity_correlation.csv"),
        ("adult SLEDAI correlation", "aSLE_severity_correlation.csv"),
        ("age-group-by-SLEDAI interaction", "age_severity_interaction.csv"),
    ]:
        table = pd.read_csv(SEVERITY / filename)
        rows.append(
            {
                "analysis": analysis,
                "interactions_tested": len(table),
                "fdr_below_0_05": int((table["fdr_q_value"] < 0.05).sum()),
                "fdr_below_0_10": int((table["fdr_q_value"] < 0.10).sum()),
            }
        )
    pd.DataFrame(rows).to_csv(MANUSCRIPT / "secondary_sledai_summary.csv", index=False)


def _style_main_table(artist, fontsize: float = 13) -> None:
    artist.auto_set_font_size(False)
    artist.set_fontsize(fontsize)
    for (row, _column), cell in artist.get_celld().items():
        cell.set_edgecolor("#B8BEC6")
        cell.set_linewidth(0.55)
        if row == 0:
            cell.set_facecolor("#DCE6F1")
            cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F5F7FA")


def render_main_tables() -> None:
    """Render the two main tables in large-type, manuscript-width layouts."""
    table1 = pd.read_csv(MANUSCRIPT / "table1_cohort_characteristics.csv")
    labels = {
        "n_donors": "Donors",
        "age_median_iqr": "Age, median (IQR)",
        "female_male": "Female/male",
        "batches_n": "Batches (n)",
        "retained_cells": "Retained cells",
        "sledai_available_n": "SLEDAI available",
        "sledai_median_iqr": "SLEDAI, median (IQR)",
    }
    transposed = table1.set_index("group").T.reset_index()
    transposed["index"] = transposed["index"].map(labels)
    transposed["index"] = transposed["index"].map(
        lambda value: _wrap_table_value(value, 16)
    )
    for column in transposed.columns[1:]:
        transposed[column] = transposed[column].map(
            lambda value: _wrap_table_value(value, 12)
        )
    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    ax.axis("off")
    artist = ax.table(
        cellText=transposed.values,
        colLabels=[
            "Characteristic",
            "Pediatric\nhealthy",
            "Pediatric\nSLE",
            "Adult\nhealthy",
            "Adult\nSLE",
        ],
        colWidths=[0.25, 0.1875, 0.1875, 0.1875, 0.1875],
        cellLoc="left",
        colLoc="left",
        loc="upper left",
        bbox=[0, 0, 1, 1],
    )
    # Match Table 2's apparent type size when both images are inserted at the
    # same manuscript width (Table 1 uses a narrower source canvas).
    _style_main_table(artist, fontsize=11)
    save_publication_figure(MANUSCRIPT / "table1_cohort_characteristics.png", fig)
    plt.close(fig)

    table2 = pd.read_csv(MANUSCRIPT / "table2_representative_interactions.csv")
    statistical = table2[
        ["Source to target", "Ligand-receptor", "Coefficient", "95% CI", "FDR"]
    ].copy()
    statistical = statistical.map(lambda value: _wrap_table_value(value, 22))
    interpretation = table2[
        [
            "Ligand-receptor",
            "Pediatric healthy to SLE",
            "Adult healthy to SLE",
            "Plain-language pattern",
        ]
    ].copy()
    interpretation = interpretation.map(lambda value: _wrap_table_value(value, 22))

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 15.5))
    for ax in axes:
        ax.axis("off")
    axes[0].set_title("A. Statistical results", fontsize=15, fontweight="bold", loc="left", pad=10)
    top = axes[0].table(
        cellText=statistical.values,
        colLabels=statistical.columns,
        colWidths=[0.28, 0.19, 0.14, 0.25, 0.14],
        cellLoc="left",
        colLoc="left",
        loc="upper left",
        bbox=[0, 0, 1, 0.94],
    )
    _style_main_table(top, fontsize=13)
    axes[1].set_title(
        "B. Group means and interpretation",
        fontsize=15,
        fontweight="bold",
        loc="left",
        pad=10,
    )
    bottom = axes[1].table(
        cellText=interpretation.values,
        colLabels=[
            "Ligand-receptor",
            "Pediatric healthy\nto SLE",
            "Adult healthy\nto SLE",
            "Plain-language pattern",
        ],
        colWidths=[0.20, 0.22, 0.22, 0.36],
        cellLoc="left",
        colLoc="left",
        loc="upper left",
        bbox=[0, 0, 1, 0.94],
    )
    _style_main_table(bottom, fontsize=13)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98, hspace=0.12)
    save_publication_figure(MANUSCRIPT / "table2_representative_interactions.png", fig)
    plt.close(fig)


def render_all_tables() -> None:
    """Create 600-DPI PNG companions for all generated manuscript CSV tables."""
    titles = {
        "table1_cohort_characteristics": "Table 1. Characteristics of the analyzed donor cohorts",
        "table2_representative_interactions": "Table 2. Representative communication patterns differing by cohort",
        "supplementary_adult_sample_exclusions": "Supplementary Table S1. Adult records excluded before analysis",
        "supplementary_cell_type_coverage": "Supplementary Table S2. Cell-type coverage by donor group",
        "supplementary_all_fdr05_interactions": "Supplementary Table S3. All interactions with FDR < 0.05",
        "significant_interaction_group_means": "Supplementary Table S4. Group means for significant interactions",
        "significant_sources": "Supplementary summary. Significant interactions by source",
        "significant_targets": "Supplementary summary. Significant interactions by target",
        "significant_source_target_pairs": "Supplementary summary. Significant source-target pairs",
        "significant_ligand_receptor_families": "Supplementary summary. Recurrent ligand-receptor pairs",
        "secondary_sledai_summary": "Supplementary summary. Secondary SLEDAI analyses",
    }
    for stem, title in titles.items():
        if stem in {"table1_cohort_characteristics", "table2_representative_interactions"}:
            continue
        table = pd.read_csv(MANUSCRIPT / f"{stem}.csv")
        save_table_png(
            table,
            MANUSCRIPT / stem,
            title=title,
            rows_per_page=16 if len(table.columns) >= 8 else 24,
            wrap_width=20 if len(table.columns) >= 8 else 26,
        )
    render_main_tables()


if __name__ == "__main__":
    cohort_table()
    cell_type_coverage()
    significant_interaction_summaries()
    adult_exclusions()
    severity_summary()
    render_all_tables()
    print(f"Wrote manuscript table CSVs and 600-DPI PNGs to {MANUSCRIPT}")
