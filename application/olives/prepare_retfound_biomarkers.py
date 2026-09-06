#!/usr/bin/env python3
"""Join central RETFound visits to the image-level OLIVES biomarker labels."""
import argparse
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, default=ROOT/'application/olives/output/retfound/manifest.csv')
    ap.add_argument('--output', type=Path, default=ROOT/'application/olives/output/retfound/visit_biomarkers.csv')
    args = ap.parse_args()
    biomarker_path = ROOT/'data/OLIVES_Dataset_Labels/ml_centric_labels/Biomarker_Clinical_Data_Images.csv'
    manifest = pd.read_csv(args.manifest, dtype={'patient_id': str, 'eye_id': str, 'visit': int})
    biomarker = pd.read_csv(biomarker_path)
    path_col = 'Path (Trial/Arm/Folder/Visit/Eye/Image Name)'
    manifest['visit_key'] = manifest['relative_path'].str.split('/').str[:5].str.join('/')
    biomarker['relative_path'] = biomarker[path_col].str.lstrip('/')
    # Biomarker labels are visit-level, but the export contains only some
    # slice filenames. Join on trial/arm/folder/visit/eye instead of basename.
    biomarker['visit_key'] = biomarker['relative_path'].str.split('/').str[:5].str.join('/')
    biomarker = biomarker.drop_duplicates('visit_key')
    selected = manifest.merge(biomarker, on='visit_key', how='left', suffixes=('', '_label'))
    label_columns = [c for c in biomarker.columns if c not in {path_col, 'relative_path', 'visit_key'}]
    missing = selected[label_columns].isna().all(axis=1)
    # The public biomarker table is sparse; preserve all manifest visits and
    # leave unavailable biomarker visits as missing for outcome-specific tests.
    # Keep one row per selected eye-visit and the 16 biomarker fields.
    keep = ['patient_id', 'eye_id', 'arm', 'study_id', 'visit', 'bcva', 'cst'] + [
        c for c in label_columns if c not in {'Eye_ID', 'Patient_ID', 'BCVA', 'CST', 'Scan (n/49)'}
    ]
    out = selected[keep].copy()
    out.to_csv(args.output, index=False)
    print(f'wrote {len(out)} visit rows with {len(keep)-7} biomarker fields; '
          f'{int((~missing).sum())} have biomarker labels')


if __name__ == '__main__':
    main()
