"""Pre-specified analysis settings shared across notebooks."""

from __future__ import annotations

RANDOM_SEED = 0
FDR_ALPHA = 0.05
EXPLORATORY_FDR_ALPHA = 0.10
MIN_CELLS_PER_SAMPLE_CELL_TYPE = 10
MIN_PATIENTS_PER_GROUP = 5
BOOTSTRAP_REPLICATES = 2000
MIN_GENES_PER_CELL = 200
MAX_MITOCHONDRIAL_PERCENT = 25.0
EXPECTED_DOUBLET_RATE = 0.06
N_HIGHLY_VARIABLE_GENES = 2000
CELLTYPIST_MAX_CELLS = 10000

PRIMARY_HYPOTHESIS = (
    "The SLE-versus-healthy difference in donor-level LIANA consensus "
    "communication scores differs between pediatric and adult participants."
)
PRIMARY_MODEL = "communication_score ~ C(disease_status) * C(age_group)"
PRIMARY_TERM = "C(disease_status)[T.SLE]:C(age_group)[T.pediatric]"

GROUP_LABELS = {
    "cH": {"age_group": "pediatric", "disease_status": "healthy"},
    "cSLE": {"age_group": "pediatric", "disease_status": "SLE"},
    "aH": {"age_group": "adult", "disease_status": "healthy"},
    "aSLE": {"age_group": "adult", "disease_status": "SLE"},
}
