#!/usr/bin/env python3
"""Build an outcome-blind TREX-DME central B-scan manifest; never alter images."""
import argparse
import collections
import hashlib
import json
import re
from pathlib import Path
from PIL import Image
from clinical_io import read_xlsx, write_csv

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path, default=ROOT/'application/olives/output/retfound')
    ap.add_argument('--max-visit', type=int, default=27)
    ap.add_argument('--min-visits', type=int, default=3)
    args = ap.parse_args()
    if args.min_visits < 1 or args.max_visit < args.min_visits:
        ap.error('Require max-visit >= min-visits >= 1')
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    workbook = ROOT/'data/OLIVES_Dataset_Labels/ml_centric_labels/Clinical_Data_Images.xlsx'
    table = read_xlsx(workbook)
    visits = collections.defaultdict(list)
    for row in table[1:]:
        d = {table[0][k]: v for k, v in row.items() if k in table[0]}
        p = Path(d.get('File_Path', ''))
        if str(p).startswith('/TREX DME/'):
            visits[(d['Eye_ID'], p.parts[-3])].append(d)
    records, exclusions = [], []
    for (eye, visit), ds in visits.items():
        vi = int(visit[1:])
        if vi > args.max_visit:
            continue
        idx = [int(re.search(r'_(\d+)\.[^.]+$', d['File_Path']).group(1)) for d in ds]
        if len(ds) != 49 or sorted(idx) != list(range(49)):
            exclusions.append(dict(eye_id=eye, visit=vi, reason='incomplete_or_ambiguous_acquisition'))
            continue
        if len({(d.get('BCVA'), d.get('CST')) for d in ds}) != 1:
            raise ValueError(f'Conflicting clinical labels: {eye}, {visit}')
        d = ds[idx.index(24)]
        p = Path(d['File_Path'])
        records.append(dict(patient_id=int(d['Patient_ID']), eye_id=int(eye), arm=p.parts[2],
                            study_id=p.parts[3], laterality=p.parts[-2], visit=vi,
                            bcva=d.get('BCVA', ''), cst=d.get('CST', ''), slice_index=24,
                            relative_path=str(p).lstrip('/')))
    by_eye = collections.defaultdict(list)
    for r in records:
        by_eye[r['eye_id']].append(r)
    candidates = collections.defaultdict(list)
    for eye, rs in by_eye.items():
        if len(rs) >= args.min_visits and any(r['visit'] == 1 for r in rs):
            candidates[rs[0]['patient_id']].append(eye)
        else:
            exclusions.append(dict(eye_id=eye, visit='all', reason='insufficient_coverage_or_missing_V1'))
    chosen = {sorted(es, key=lambda e: (-len(by_eye[e]), e))[0] for es in candidates.values()}
    for es in candidates.values():
        for eye in set(es)-chosen:
            exclusions.append(dict(eye_id=eye, visit='all', reason='fellow_study_eye'))
    records = sorted((r for r in records if r['eye_id'] in chosen), key=lambda r:(r['patient_id'], r['visit']))
    seen, retained = {}, []
    for r in records:
        path = ROOT/'data/OLIVES'/r['relative_path']
        with Image.open(path) as im:
            im = im.convert('L')
            r['height'], r['width'] = im.height, im.width
            digest = hashlib.sha256(im.tobytes()).hexdigest()
        key = (r['eye_id'], r['height'], r['width'], digest)
        if key in seen:
            exclusions.append(dict(eye_id=r['eye_id'], visit=r['visit'], reason=f'exact_pixel_copy_of_V{seen[key]}'))
            continue
        seen[key] = r['visit']
        r['pixel_sha256'] = digest
        r['feature_row'] = len(retained)
        retained.append(r)
    write_csv(out/'manifest.csv', retained)
    write_csv(out/'exclusions.csv', exclusions)
    summary = dict(n_patients=len(candidates), n_images=len(retained), cohort='TREX-DME only',
                   time_axis='documented visit index, not elapsed weeks', max_visit=args.max_visit,
                   min_visits=args.min_visits, selection='V1 required; one eye per patient by coverage then Eye_ID',
                   image_processing='none; original central B-scan passed to official model evaluation transform',
                   clinical_workbook_sha256=hashlib.sha256(workbook.read_bytes()).hexdigest())
    (out/'cohort.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
