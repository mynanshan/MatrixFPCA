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

TRAIN subjects are ordered by a stratified random ranking (`train_rank`) and
the first --n-train are exported. The strata cross EF quintile with eight
image-appearance clusters (k-means on PCA scores of each video's
time-averaged image and ES - ED difference image, 28 x 28), so every kind of
heart and every kind of image is represented in proportion, including the
rare low-EF hearts. Because the ranking is systematic within strata, any
prefix (first 250, 500, 1000, ...) is itself a proportional stratified
sample, and ranks 501-1000 form a second, disjoint one. The ranking depends
only on the preprocessing output and --seed, not on --pool.

The representativeness of the exported TRAIN sample against all retained
TRAIN videos is written to sample_check.csv, sample_strata.csv and
sample_representativeness.png. VAL and TEST are random samples (or all).
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


def appearance_features(prep, v, mask, seed, n_pc=20):
    """PCA scores of each TRAIN video's time-averaged and ES - ED images
    (gain-normalised, 4 x 4 pooled, masked). Cached in the prep directory."""
    cache = os.path.join(prep, "appearance_features.npz")
    ids = v.id.values
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        if list(z["ids"]) == list(ids):
            return z["scores"]
    m4 = pool_frames(mask[None].astype(np.float32), 4)[0] >= 1.0
    F = np.empty((len(ids), 2 * int(m4.sum())), np.float32)
    for i, r in enumerate(v.itertuples()):
        f = np.load(os.path.join(prep, "frames", r.id + ".npz"))["frames"]
        f = pool_frames(f.astype(np.float32) / 255.0 / r.norm_scale, 4)
        F[i] = np.r_[f.mean(0)[m4], (f[-1] - f[0])[m4]]
    F -= F.mean(0)
    u, s, _ = np.linalg.svd(F, full_matrices=False)
    scores = u[:, :n_pc] * s[:n_pc]
    np.savez(cache, ids=ids, scores=scores)
    return scores


def stratified_rank(d, scores, seed, n_clusters=8, n_ef=5):
    """Systematic stratified random ranking of TRAIN videos.

    Strata = EF quintile x appearance cluster. Within a stratum of size N_s
    the videos get random positions (i + U) / N_s, i = 0..N_s-1; sorting all
    videos by position makes every prefix a proportional stratified sample."""
    from scipy.cluster.vq import kmeans2
    rng = np.random.default_rng(seed)
    z = scores[:, :10] / scores[:, :10].std(0)
    _, cl = kmeans2(z, n_clusters, seed=seed, minit="++")
    ef_q = pd.qcut(d.EF.values, n_ef, labels=False)
    stratum = ef_q * n_clusters + cl
    pos = np.empty(len(d))
    for s in np.unique(stratum):
        idx = np.flatnonzero(stratum == s)
        pos[idx] = (rng.permutation(len(idx)) + rng.random(len(idx))) / len(idx)
    order = np.lexsort((rng.random(len(d)), pos))
    rank = np.empty(len(d), int)
    rank[order] = np.arange(1, len(d) + 1)
    return rank, cl, ef_q


def ks_stat(a, b):
    a, b = np.sort(a), np.sort(b)
    x = np.r_[a, b]
    return np.abs(np.searchsorted(a, x, "right") / len(a)
                  - np.searchsorted(b, x, "right") / len(b)).max()


SAMPLE_VARS = {
    "EF": "ejection fraction (%)",
    "EDV": "end-diastolic volume (mL)",
    "ESV": "end-systolic volume (mL)",
    "norm_scale": "image gain (mean in-mask grey)",
    "axis_len": "native LV length (px, before alignment)",
    "axis_tilt": "native LV tilt (deg, before alignment)",
    "FPS": "frame rate (fps)",
    "gap_s": "systole duration ED->ES (s)",
}


def sample_check(pop, smp, scores, out):
    """Compare the TRAIN sample with all retained TRAIN videos."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for c, label in SAMPLE_VARS.items():
        a, b = pop[c].values, smp[c].values
        q = lambda x: np.percentile(x, [1, 5, 50, 95, 99])
        rows.append(dict(variable=c, description=label,
                         **{f"all_p{k}": v for k, v in zip((1, 5, 50, 95, 99), q(a))},
                         **{f"sample_p{k}": v for k, v in zip((1, 5, 50, 95, 99), q(b))},
                         ks_D=ks_stat(a, b),
                         ks_D_crit_5pct=1.36 * np.sqrt(1 / len(a) + 1 / len(b))))
    for lo, hi, label in ((0, 30, "EF < 30"), (0, 40, "EF < 40"), (40, 55, "EF 40-55"),
                          (75, 101, "EF >= 75")):
        e = ((pop.EF >= lo) & (pop.EF < hi)).mean() * len(smp)
        g = ((smp.EF >= lo) & (smp.EF < hi)).sum()
        rows.append(dict(variable=label, description="subjects: expected vs sampled",
                         all_p50=e, sample_p50=g))
    pd.DataFrame(rows).to_csv(os.path.join(out, "sample_check.csv"), index=False)
    st = (pop.groupby(["ef_quintile", "appearance_cluster"]).size().rename("all").to_frame()
          .join(smp.groupby(["ef_quintile", "appearance_cluster"]).size().rename("sample"))
          .fillna(0).reset_index())
    st["expected"] = st["all"] * len(smp) / len(pop)
    st.to_csv(os.path.join(out, "sample_strata.csv"), index=False)

    fig, ax = plt.subplots(2, 5, figsize=(17, 6.4))
    for a_, c in zip(ax.flat[:8], SAMPLE_VARS):
        for d, col, lab in ((pop, "#9aa5b1", f"all TRAIN (n={len(pop)})"),
                            (smp, "#1f3f8f", f"sample (n={len(smp)})")):
            x = np.sort(d[c].values)
            a_.step(x, np.arange(1, len(x) + 1) / len(x), color=col, lw=1.6, label=lab)
        a_.set_title(f"{SAMPLE_VARS[c]}\nKS D = {ks_stat(pop[c].values, smp[c].values):.3f}",
                     fontsize=9)
        a_.set_ylabel("cumulative share", fontsize=8)
    ax.flat[0].legend(fontsize=7, frameon=False, loc="lower right")
    inp = pop.index.isin(smp.index)
    a_ = ax.flat[8]
    a_.scatter(scores[~inp, 0], scores[~inp, 1], s=3, c="#c8ccd2", label="not sampled")
    a_.scatter(scores[inp, 0], scores[inp, 1], s=5, c="#1f3f8f", label="sample")
    a_.set_title("image appearance (PCA of time-average\nand ES - ED images)", fontsize=9)
    a_.set_xlabel("appearance PC1", fontsize=8); a_.set_ylabel("appearance PC2", fontsize=8)
    a_.legend(fontsize=7, frameon=False, markerscale=3)
    a_ = ax.flat[9]
    cl = st.groupby("appearance_cluster")[["all", "sample"]].sum()
    x = np.arange(len(cl))
    a_.bar(x - 0.2, cl["all"] / cl["all"].sum(), 0.4, color="#9aa5b1", label="all TRAIN")
    a_.bar(x + 0.2, cl["sample"] / cl["sample"].sum(), 0.4, color="#1f3f8f", label="sample")
    a_.set_xticks(x, [str(i + 1) for i in cl.index])
    a_.set_xlabel("appearance cluster", fontsize=8); a_.set_ylabel("share", fontsize=8)
    a_.set_title("appearance clusters", fontsize=9)
    a_.legend(fontsize=7, frameon=False)
    for a_ in ax.flat:
        a_.spines[["top", "right"]].set_visible(False)
        a_.tick_params(labelsize=7)
    fig.suptitle(f"Is the TRAIN sample representative? (n = {len(smp)} of {len(pop)})")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "sample_representativeness.png"), dpi=90)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--prep", default="application/echonet/output/prep_similarity")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", type=int, default=2)
    ap.add_argument("--n-train", type=int, default=500, help="0 = all")
    ap.add_argument("--n-val", type=int, default=300, help="0 = all")
    ap.add_argument("--n-test", type=int, default=0, help="0 = all")
    ap.add_argument("--normalise", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20261124)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    v = pd.read_csv(os.path.join(a.prep, "videos.csv"))
    v = v[v.exclusion.isna() | (v.exclusion == "")]
    mask = np.load(os.path.join(a.prep, "mask.npy"))
    tr = v[v.Split == "TRAIN"].sort_values("id").reset_index(drop=True)
    scores = appearance_features(a.prep, tr, mask, a.seed)
    tr["train_rank"], cl, efq = stratified_rank(tr, scores, a.seed)
    tr["appearance_cluster"], tr["ef_quintile"] = cl, efq
    n_tr = len(tr) if a.n_train <= 0 else min(a.n_train, len(tr))
    smp = tr[tr.train_rank <= n_tr]
    sample_check(tr, smp, scores, a.out)
    parts = [smp]
    for split, n in (("VAL", a.n_val), ("TEST", a.n_test)):
        d = v[v.Split == split]
        parts.append(d if n <= 0 or n >= len(d) else d.sample(n, random_state=a.seed))
    s = pd.concat(parts).sort_values(["Split", "id"]).reset_index(drop=True)

    mask_p = pool_frames(mask[None].astype(np.float32), a.pool)[0] >= 1.0   # all sub-pixels in
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
         "norm_scale", "colour_frac", "axis_len", "axis_tilt", "n_frames",
         "train_rank", "appearance_cluster", "ef_quintile"]
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
