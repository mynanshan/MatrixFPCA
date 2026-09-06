#!/usr/bin/env python3
"""Extract reusable patient-level clinical quantities from the OLIVES DME workbook."""
import argparse
import csv
from pathlib import Path
from clinical_io import read_xlsx, write_csv

ROOT = Path(__file__).resolve().parents[2]
WANTED = [
    'Age at Enrollment', 'Gender', 'Ethnicity', 'Height (In.)', 'Height (in.)',
    'Weight', 'Weight (kg)', 'BMI', 'Dia-Type', 'Dia-Duration', 'InsulinUser',
    'Smoker', 'SmokingDuration', 'CataractSurgery', 'Bi-Lateral', 'DR Level',
    '# of Tx visits in first year', '# of Tx visits in 24 months',
    'Change in ETDRS', 'Change in ETDRS @ 2 years',
    '1-year 1 line gainers (0=No, 1= Yes)', '1-year 1 line losers (0=No, 1= Yes)',
    '1-year 3 line gainers (0=No, 1 = Yes)', '1-year 3 line losers (0=No, 1 = Yes)',
]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, default=ROOT/'application/olives/output/retfound/manifest.csv')
    ap.add_argument('--output', type=Path, default=ROOT/'application/olives/output/retfound/demographics.csv')
    args = ap.parse_args()
    workbook = ROOT/'data/OLIVES_Dataset_Labels/full_labels/OCT-DME.xlsx'
    patient_values = {}
    for sheet in (1, 2, 3):
        rows = read_xlsx(workbook, sheet)
        header = rows[3] if len(rows) > 3 else {}
        for row in rows[4:]:
            study = row.get('FD') or row.get('A')
            if not study:
                continue
            values = {name: row.get(key, '') for key, name in header.items() if name in WANTED}
            values.setdefault('Age at Enrollment', row.get('FE', ''))
            patient_values[str(study)] = values
    with args.manifest.open(newline='') as f:
        manifest = list(csv.DictReader(f))
    records = []
    seen = set()
    for row in manifest:
        patient, study = row['patient_id'], row['study_id']
        if patient in seen:
            continue
        seen.add(patient)
        values = patient_values.get(study, {})
        record = {'patient_id': patient, 'study_id': study}
        record['age_at_enrollment'] = values.get('Age at Enrollment', '')
        for name in WANTED:
            if name == 'Age at Enrollment':
                continue
            record[name] = values.get(name, '')
        records.append(record)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_csv(args.output, records)
    print(f'wrote {len(records)} patient demographic/clinical rows to {args.output}')

if __name__ == '__main__':
    main()

