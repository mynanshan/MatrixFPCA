#!/usr/bin/env python3
"""Screen baseline RETFound embeddings against a broad clinical panel.

The PCA is fitted only to baseline frozen embeddings.  Baseline scores are
then tested against baseline clinical quantities and V13-minus-V1 changes;
the clinical labels do not enter the PCA fit.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA


EXCLUDE_BIOMARKER = {"patient_id", "eye_id", "arm", "study_id", "visit", "bcva", "cst"}
EXCLUDE_DEMOGRAPHIC = {
    "patient_id", "study_id", "age_at_enrollment", "Gender", "Ethnicity",
    "Dia-Type", "DR Level",
}


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def bh(values):
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    positions = np.flatnonzero(finite)
    order = positions[np.argsort(values[finite])]
    running = 1.0
    for rank in range(len(order) - 1, -1, -1):
        position = order[rank]
        running = min(running, values[position] * len(order) / (rank + 1))
        out[position] = running
    return out


def adjusted_corr(x, y, arm):
    keep = np.isfinite(x) & np.isfinite(y)
    x, y, arm = np.asarray(x)[keep], np.asarray(y)[keep], np.asarray(arm)[keep]
    if len(x) < 8 or len(np.unique(arm)) < 2:
        return np.nan, np.nan, len(x)
    design = np.column_stack([np.ones(len(arm))] + [arm == a for a in sorted(set(arm))[1:]])
    xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(xr) == 0 or np.std(yr) == 0:
        return np.nan, np.nan, len(x)
    r, p = pearsonr(xr, yr)
    return float(r), float(p), len(x)


def outcome_screenable(values, arms):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    vals = values[finite]
    if len(vals) < 10 or len(np.unique(arms[finite])) < 2 or len(np.unique(vals)) < 2:
        return False, len(vals), len(np.unique(vals)), 0
    unique, counts = np.unique(vals, return_counts=True)
    min_level = int(counts.min())
    # Very low-cardinality labels need at least five observations in each level.
    return (len(unique) > 5 or min_level >= 5), len(vals), len(unique), min_level


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=root / "application/olives/output/retfound/features")
    parser.add_argument("--representation", choices=("cls", "pooled"), default="cls")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--demographics", type=Path)
    parser.add_argument("--biomarkers", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--components", type=int, default=30)
    args = parser.parse_args()

    manifest_path = args.manifest or args.input_dir.parent / "manifest.csv"
    demographics_path = args.demographics or args.input_dir.parent / "demographics.csv"
    biomarkers_path = args.biomarkers or args.input_dir.parent / "visit_biomarkers.csv"
    output_dir = args.output_dir or args.input_dir.parent / "embedding_screen"
    output_dir.mkdir(parents=True, exist_ok=True)

    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle))
    embeddings = np.load(args.input_dir / f"{args.representation}.npy", mmap_mode="r")
    if embeddings.ndim != 2 or embeddings.shape[0] != len(manifest):
        raise ValueError(f"Expected observation-by-feature matrix, got {embeddings.shape}")

    patients = sorted({row["patient_id"] for row in manifest})
    by_patient = {p: [i for i, row in enumerate(manifest) if row["patient_id"] == p] for p in patients}
    baseline_index = []
    v13_index = {}
    arms = []
    for patient in patients:
        indices = by_patient[patient]
        baseline_index.append(next(i for i in indices if int(manifest[i]["visit"]) == 1))
        target = next((i for i in indices if int(manifest[i]["visit"]) == 13), None)
        if target is not None:
            v13_index[patient] = target
        arms.append(manifest[baseline_index[-1]]["arm"])
    baseline_index = np.asarray(baseline_index, dtype=int)
    arms = np.asarray(arms)
    n_components = min(args.components, len(patients) - 1, embeddings.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full")
    pca.fit(np.asarray(embeddings[baseline_index], dtype=np.float64))
    baseline_scores = pca.transform(np.asarray(embeddings[baseline_index], dtype=np.float64))

    with demographics_path.open(newline="") as handle:
        demographics_rows = {row["patient_id"]: row for row in csv.DictReader(handle)}
    with biomarkers_path.open(newline="") as handle:
        biomarker_rows = [row for row in csv.DictReader(handle) if row["visit"] in {"1", "1.0"}]
    biomarkers = {row["patient_id"]: row for row in biomarker_rows}
    biomarker_names = sorted(set().union(*(row.keys() for row in biomarkers.values())) - EXCLUDE_BIOMARKER)

    outcomes = {
        "age": [numeric(demographics_rows.get(p, {}).get("age_at_enrollment")) for p in patients],
        "baseline_bcva": [numeric(manifest[i]["bcva"]) for i in baseline_index],
        "baseline_cst": [numeric(manifest[i]["cst"]) for i in baseline_index],
        "v13_bcva_change": [],
        "v13_cst_change": [],
    }
    for patient, i0 in zip(patients, baseline_index):
        i13 = v13_index.get(patient)
        outcomes["v13_bcva_change"].append(
            numeric(manifest[i13]["bcva"]) - numeric(manifest[i0]["bcva"]) if i13 is not None else np.nan
        )
        outcomes["v13_cst_change"].append(
            numeric(manifest[i13]["cst"]) - numeric(manifest[i0]["cst"]) if i13 is not None else np.nan
        )
    for name in biomarker_names:
        outcomes[f"baseline biomarker: {name}"] = [numeric(biomarkers.get(p, {}).get(name)) for p in patients]

    demographic_names = sorted(set().union(*(row.keys() for row in demographics_rows.values())) - EXCLUDE_DEMOGRAPHIC)
    for name in demographic_names:
        values = [numeric(demographics_rows.get(p, {}).get(name)) for p in patients]
        finite = np.asarray(values)[np.isfinite(values)]
        if len(finite) >= 10 and len(np.unique(finite)) > 1:
            outcomes[f"demographic: {name}"] = values

    records = []
    for component in range(n_components):
        for outcome, values in outcomes.items():
            y = np.asarray(values, dtype=float)
            r, p, n = adjusted_corr(baseline_scores[:, component], y, arms)
            keep = np.isfinite(baseline_scores[:, component]) & np.isfinite(y)
            if keep.sum() >= 5 and np.std(y[keep]) > 0:
                rho, rho_p = spearmanr(baseline_scores[keep, component], y[keep])
                rho, rho_p = float(rho), float(rho_p)
            else:
                rho, rho_p = np.nan, np.nan
            screenable, outcome_n, unique_n, min_level_n = outcome_screenable(y, arms)
            records.append({
                "mode": component + 1,
                "outcome": outcome,
                "n": n,
                "pve": float(pca.explained_variance_ratio_[component]),
                "arm_adjusted_r": r,
                "arm_adjusted_p": p,
                "spearman_rho": rho,
                "spearman_p": rho_p,
                "outcome_n": outcome_n,
                "outcome_unique": unique_n,
                "outcome_min_level_n": min_level_n,
                "outcome_screenable": screenable,
            })

    p = np.asarray([row["arm_adjusted_p"] for row in records], dtype=float)
    rho_p = np.asarray([row["spearman_p"] for row in records], dtype=float)
    q = bh(p)
    rho_q = bh(rho_p)
    screenable_p = np.where([row["outcome_screenable"] for row in records], p, np.nan)
    screenable_q = bh(screenable_p)
    for row, q_value, rho_value, screen_q in zip(records, q, rho_q, screenable_q):
        row["arm_adjusted_q"] = float(q_value) if np.isfinite(q_value) else np.nan
        row["spearman_q"] = float(rho_value) if np.isfinite(rho_value) else np.nan
        row["arm_adjusted_q_screenable_global"] = float(screen_q) if np.isfinite(screen_q) else np.nan
    for outcome in outcomes:
        indices = [i for i, row in enumerate(records) if row["outcome"] == outcome]
        outcome_q = bh(p[indices])
        for i, value in zip(indices, outcome_q):
            records[i]["arm_adjusted_q_within_outcome"] = float(value) if np.isfinite(value) else np.nan

    fields = list(records[0])
    with (output_dir / f"{args.representation}_baseline_pca_clinical_associations.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    with (output_dir / f"{args.representation}_baseline_pca_spectrum.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mode", "pve", "cumulative_pve"])
        writer.writeheader()
        cumulative = 0.0
        for mode, pve in enumerate(pca.explained_variance_ratio_, 1):
            cumulative += float(pve)
            writer.writerow({"mode": mode, "pve": float(pve), "cumulative_pve": cumulative})
    score_fields = ["patient_id", "arm"] + [f"mode_{j}" for j in range(1, n_components + 1)]
    with (output_dir / f"{args.representation}_baseline_pca_scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=score_fields)
        writer.writeheader()
        for i, patient in enumerate(patients):
            row = {"patient_id": patient, "arm": arms[i]}
            row.update({f"mode_{j}": float(baseline_scores[i, j - 1]) for j in range(1, n_components + 1)})
            writer.writerow(row)
    report = {
        "representation": args.representation,
        "n_patients": len(patients),
        "n_baseline": len(baseline_index),
        "n_components": n_components,
        "n_outcomes": len(outcomes),
        "n_tests": len(records),
        "screenable_outcomes": len({r["outcome"] for r in records if r["outcome_screenable"]}),
        "screenable_tests": sum(r["outcome_screenable"] for r in records),
        "interpretation": "Unsupervised baseline PCA; clinical labels are used only after the PCA fit.",
    }
    (output_dir / f"{args.representation}_baseline_pca_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
