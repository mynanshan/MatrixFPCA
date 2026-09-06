#!/usr/bin/env python3
"""Screen frozen RETFound visit features for patient-level clinical signal."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    out = np.empty_like(values)
    running = 1.0
    for pos in range(len(values) - 1, -1, -1):
        rank = pos + 1
        running = min(running, values[order[pos]] * len(values) / rank)
        out[order[pos]] = running
    return out


def adjusted_corr(x, y, arm):
    keep = np.isfinite(x) & np.isfinite(y)
    x, y, arm = x[keep], y[keep], arm[keep]
    if len(x) < 5 or len(np.unique(arm)) < 1:
        return np.nan, np.nan, len(x)
    design = np.column_stack([np.ones(len(x))] + [arm == a for a in np.unique(arm)[1:]])
    xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    r, p = pearsonr(xr, yr)
    return float(r), float(p), len(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir', type=Path, default=Path('application/olives/output/retfound/features'))
    ap.add_argument('--manifest', type=Path)
    ap.add_argument('--components', type=int, default=10)
    args = ap.parse_args()
    manifest = args.manifest or args.input_dir.parent / 'manifest.csv'
    with manifest.open(newline='') as f:
        rows = list(csv.DictReader(f))
    pooled = np.load(args.input_dir / 'pooled.npy', mmap_mode='r')
    if pooled.shape != (len(rows), 1024):
        raise ValueError(f'Expected pooled shape {(len(rows), 1024)}, got {pooled.shape}')
    patients = sorted({r['patient_id'] for r in rows})
    by_patient = {}
    for patient in patients:
        ix = [i for i, r in enumerate(rows) if r['patient_id'] == patient]
        by_patient[patient] = ix
    # The PCA is an unsupervised summary of the frozen representation. It is
    # used only to make a low-dimensional association table, not to fit outcomes.
    baseline_ix = [next(i for i in ix if int(rows[i]['visit']) == 1) for ix in by_patient.values()]
    pca = PCA(n_components=min(args.components, len(baseline_ix) - 1), svd_solver='full')
    pca.fit(np.asarray(pooled[baseline_ix], dtype=np.float64))
    scores = pca.transform(np.asarray(pooled, dtype=np.float64))
    records = []
    outcomes = []
    base = {patient: i for patient, i in zip(patients, baseline_ix)}
    for name, visit, delta in [('baseline_bcva', 1, False), ('baseline_cst', 1, False),
                               ('v13_bcva_change', 13, True), ('v13_cst_change', 13, True)]:
        vals = []
        feature_rows = []
        arms = []
        used = []
        for p, ix in by_patient.items():
            i0 = base[p]
            target = next((i for i in ix if int(rows[i]['visit']) == visit), None)
            if target is None:
                continue
            field = 'bcva' if 'bcva' in name else 'cst'
            try:
                y0, y1 = float(rows[i0][field]), float(rows[target][field])
            except (TypeError, ValueError):
                continue
            if not np.isfinite(y0) or not np.isfinite(y1):
                continue
            vals.append(y1 - y0 if delta else y0)
            feature_rows.append(scores[target] - scores[i0] if delta else scores[i0])
            arms.append(rows[i0]['arm'])
            used.append(p)
        if len(vals) < 5:
            continue
        outcomes.append(name)
        y = np.asarray(vals)
        X = np.asarray(feature_rows)
        arm = np.asarray(arms)
        for component in range(X.shape[1]):
            rho, rho_p = spearmanr(X[:, component], y)
            r_adj, p_adj, n = adjusted_corr(X[:, component], y, arm)
            records.append(dict(outcome=name, component=component + 1, n=n,
                                spearman_rho=float(rho), spearman_p=float(rho_p),
                                arm_adjusted_r=float(r_adj), arm_adjusted_p=float(p_adj),
                                explained_variance=float(pca.explained_variance_ratio_[component])))
    for field in ('spearman_p', 'arm_adjusted_p'):
        q = bh([r[field] for r in records])
        for r, value in zip(records, q):
            r[field.replace('_p', '_q')] = float(value)
    args.input_dir.mkdir(parents=True, exist_ok=True)
    with (args.input_dir / 'clinical_feature_associations.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]) if records else ['outcome'])
        writer.writeheader(); writer.writerows(records)
    report = dict(n_patients=len(patients), n_baseline=len(baseline_ix), outcomes=outcomes,
                  pca_components=pca.n_components_, explained_variance=pca.explained_variance_ratio_.tolist(),
                  interpretation='Exploratory patient-level association screen; PCA and tests are descriptive.')
    (args.input_dir / 'clinical_feature_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
