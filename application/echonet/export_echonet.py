#!/usr/bin/env python
"""Export preprocessed EchoNet frames for R (step 1 of the RDS conversion).

Selects subjects, block-averages frames (--pool), divides each video by its
gain-normalisation scale, applies the analysis mask, and writes

  frames.f32     float32, one frame after another, each frame column-major
                 (so R's matrix(x, p, q) restores the image)
  frames.csv     one row per frame: id, split, grid_idx (0-based), phase,
                 frame_idx, offset (0-based frame position in frames.f32)
  subjects.csv   one row per subject: id, split, EF, ESV, EDV, n_frames, ...
  mask.csv       p x q 0/1 analysis mask at the pooled resolution
  export.json    settings

`echonet_to_rds.R` then turns this into one rTensor RDS per subject.

TRAIN subjects are sampled stratified by EF decile, so the rare low-EF hearts
are represented in proportion; VAL and TEST are random samples (or all).
"""
import argparse
import json
import os

import numpy as np
import pandas as pd


def pool_frames(f, k):
    if k == 1:
        return f
    m, p, q = f.shape
    return f.reshape(m, p // k, k, q // k, k).mean(axis=(2, 4))


def stratified_sample(d, n, seed):
    if n <= 0 or n >= len(d):
        return d
    dec = pd.qcut(d.EF, 10, labels=False, duplicates="drop")
    per = np.floor(n * dec.value_counts(normalize=True)).astype(int)
    out = [d[dec == g].sample(per[g], random_state=seed) for g in per.index]
    out = pd.concat(out)
    rest = d.drop(out.index)
    if len(out) < n:
        out = pd.concat([out, rest.sample(n - len(out), random_state=seed)])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--prep", default="application/echonet/output/prep_similarity")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", type=int, default=2)
    ap.add_argument("--n-train", type=int, default=500)
    ap.add_argument("--n-val", type=int, default=300, help="0 = all")
    ap.add_argument("--n-test", type=int, default=0, help="0 = all")
    ap.add_argument("--normalise", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20261124)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    v = pd.read_csv(os.path.join(a.prep, "videos.csv"))
    v = v[v.exclusion.isna() | (v.exclusion == "")]
    parts = []
    for split, n in (("TRAIN", a.n_train), ("VAL", a.n_val), ("TEST", a.n_test)):
        d = v[v.Split == split]
        parts.append(stratified_sample(d, n, a.seed) if split == "TRAIN"
                     else (d if n <= 0 or n >= len(d)
                           else d.sample(n, random_state=a.seed)))
    s = pd.concat(parts).sort_values(["Split", "id"]).reset_index(drop=True)

    mask = np.load(os.path.join(a.prep, "mask.npy")).astype(np.float32)
    mask_p = pool_frames(mask[None], a.pool)[0] >= 1.0   # all sub-pixels in
    p, q = mask_p.shape

    rows, off = [], 0
    with open(os.path.join(a.out, "frames.f32"), "wb") as fh:
        for r in s.itertuples():
            z = np.load(os.path.join(a.prep, "frames", r.id + ".npz"))
            f = z["frames"].astype(np.float32) / 255.0
            if a.normalise:
                f /= r.norm_scale
            f = pool_frames(f, a.pool) * mask_p
            for j in range(len(f)):
                fh.write(np.ascontiguousarray(f[j].T).tobytes())
            rows.append(pd.DataFrame(dict(
                id=r.id, split=r.Split, grid_idx=z["grid_idx"],
                phase=z["phase"], frame_idx=z["frame_idx"],
                offset=np.arange(off, off + len(f)))))
            off += len(f)
    fr = pd.concat(rows, ignore_index=True)
    fr.to_csv(os.path.join(a.out, "frames.csv"), index=False)
    s["n_frames"] = s.id.map(fr.groupby("id").size())
    s.rename(columns={"Split": "split"})[
        ["id", "split", "EF", "ESV", "EDV", "FPS", "ED", "ES", "gap_s",
         "norm_scale", "colour_frac", "axis_len", "axis_tilt", "n_frames"]
    ].to_csv(os.path.join(a.out, "subjects.csv"), index=False)
    np.savetxt(os.path.join(a.out, "mask.csv"), mask_p.astype(int), fmt="%d",
               delimiter=",")
    with open(os.path.join(a.out, "export.json"), "w") as fh:
        json.dump(dict(vars(a), p=p, q=q, n_frames=off,
                       mask_pixels=int(mask_p.sum()),
                       n_subjects=s.Split.value_counts().to_dict()), fh, indent=2)
    print(f"exported {len(s)} subjects, {off} frames, {p}x{q}, "
          f"mask {mask_p.sum()} px -> {a.out}")
    print(s.Split.value_counts().to_string())


if __name__ == "__main__":
    main()
