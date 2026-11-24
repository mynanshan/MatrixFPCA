#!/usr/bin/env python
"""Go/no-go exploration of preprocessed EchoNet data (PLAN.md Section 4).

Reads one or more preprocessing directories (one per alignment mode) and
reports, for each:
  * variance decomposition of in-mask pixels on the registered systolic phase
    grid: between-subject vs within-subject (temporal), and the part of the
    temporal variance explained by the mean dynamics;
  * EF prediction from simple (non-functional) PCA baselines: PCA on a
    TRAIN sample, ridge regression with the penalty tuned on VAL, evaluated
    on TEST (R^2, MAE, AUC for EF < 40%).

Frames are linearly interpolated in phase onto an 11-point grid
(t = 0, 0.1, ..., 1) so that every subject has the same design here. This is
only for the cheap checks; MatFPCA itself uses the irregular frames.

Writes tables (csv + markdown) and figures to --report.
"""
import argparse
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PH = np.linspace(0, 1, 11)


def load_subject(prep, vid):
    z = np.load(os.path.join(prep, "frames", vid + ".npz"))
    f = z["frames"].astype(np.float32) / 255.0
    ph = z["phase"]
    idx = np.searchsorted(ph, PH).clip(1, len(ph) - 1)
    w = ((PH - ph[idx - 1]) / (ph[idx] - ph[idx - 1]))[:, None, None]
    return (1 - w) * f[idx - 1] + w * f[idx]


def load_split(prep, vids, mask):
    return np.stack([load_subject(prep, v)[:, mask] for v in vids])


def var_decomp(X):
    """X: n x T x P. Shares of total variance around the grand mean."""
    gm = X.mean(axis=(0, 1))
    subj = X.mean(1)
    between = ((subj - gm) ** 2).mean()
    within = ((X - subj[:, None]) ** 2).mean()
    common = ((X.mean(0) - gm) ** 2).mean()   # mean dynamics mu(t) - mu
    tot = between + within
    return dict(between=between / tot, within=within / tot,
                mean_dynamics=common / tot,
                subject_dynamics=(within - common) / tot)


def auc(y, s):
    pos, neg = s[y], s[~y]
    r = pd.Series(np.r_[pos, neg]).rank().values
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def ridge_eval(Ftr, ytr, Fva, yva, Fte, yte, k):
    mu = Ftr.mean(0)
    _, _, Vt = np.linalg.svd(Ftr - mu, full_matrices=False)
    P = Vt[:k].T
    Z = [(F - mu) @ P for F in (Ftr, Fva, Fte)]
    sd = Z[0].std(0)
    Z = [z / sd for z in Z]
    best = None
    for lam in 10.0 ** np.arange(-2, 4.5, 0.5):
        A = np.c_[np.ones(len(ytr)), Z[0]]
        pen = np.diag([0] + [lam] * k)
        b = np.linalg.solve(A.T @ A + pen, A.T @ ytr)
        rmse = np.sqrt(np.mean((np.c_[np.ones(len(yva)), Z[1]] @ b - yva) ** 2))
        if best is None or rmse < best[0]:
            best = (rmse, lam, b)
    b = best[2]
    pred = np.c_[np.ones(len(yte)), Z[2]] @ b
    r2 = 1 - ((yte - pred) ** 2).sum() / ((yte - yte.mean()) ** 2).sum()
    return dict(R2=r2, MAE=np.abs(yte - pred).mean(),
                AUC_EF40=auc(yte < 40, -pred), lambda_=best[1])


def features(X):
    n = len(X)
    return {
        "ED frame (static)": X[:, 0],
        "Time average": X.mean(1),
        "ES - ED difference": X[:, -1] - X[:, 0],
        "Stacked video": X.reshape(n, -1),
    }


def maps_figure(X, mask, title, path):
    def img(v):
        out = np.full(mask.shape, np.nan)
        out[mask] = v
        return out
    mu = X.mean(0)
    fig, ax = plt.subplots(2, 5, figsize=(14, 6.2))
    for k, j in enumerate([0, 3, 5, 8, 10]):
        ax[0, k].imshow(img(mu[j]), cmap="gray", vmin=0, vmax=0.6)
        ax[0, k].set_title(f"mean frame, phase {PH[j]:.1f}")
    d = mu[-1] - mu[0]
    v = np.nanmax(np.abs(d))
    ax[1, 0].imshow(img(d), cmap="RdBu_r", vmin=-v, vmax=v)
    ax[1, 0].set_title("mean ES − mean ED")
    sdb = X.mean(1).std(0)
    sdt = X.std(1).mean(0)
    vmax = max(sdb.max(), sdt.max())
    ax[1, 1].imshow(img(sdb), cmap="magma", vmin=0, vmax=vmax)
    ax[1, 1].set_title("between-subject SD")
    ax[1, 2].imshow(img(sdt), cmap="magma", vmin=0, vmax=vmax)
    ax[1, 2].set_title("within-subject temporal SD")
    ax[1, 3].imshow(img(sdt / sdb), cmap="viridis", vmin=0, vmax=1.2)
    ax[1, 3].set_title("temporal / between SD")
    ax[1, 4].imshow(mask, cmap="gray")
    ax[1, 4].set_title(f"analysis mask ({mask.sum()} px)")
    for a in ax.flat:
        a.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def ef_figure(et, vt, path):
    """Bar charts: EF test R^2 by feature set and alignment (k = 50), and the
    variance shares by alignment."""
    modes = list(dict.fromkeys(et.alignment))
    feats = list(dict.fromkeys(et.features))
    colours = ["#9aa5b1", "#5b8def", "#1f3f8f"][:len(modes)]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2),
                           gridspec_kw=dict(width_ratios=[1.7, 1]))
    kmax = et.k.max()
    w = 0.8 / len(modes)
    for j, m in enumerate(modes):
        d = et[(et.alignment == m) & (et.k == kmax)].set_index("features").loc[feats]
        x = np.arange(len(feats)) + (j - (len(modes) - 1) / 2) * w
        ax[0].bar(x, d.R2, w, color=colours[j], label=f"alignment: {m}")
        for xi, r in zip(x, d.R2):
            ax[0].text(xi, r + 0.005, f"{r:.2f}", ha="center", fontsize=8)
    ax[0].set_xticks(np.arange(len(feats)), feats)
    ax[0].set_ylabel(f"EF test R² (ridge on {kmax} PCs)")
    ax[0].set_title("EF is carried by the motion, not the static frame")
    ax[0].legend(frameon=False, fontsize=9)
    ax[0].spines[["top", "right"]].set_visible(False)
    d = vt[~vt.normalised].set_index("alignment").loc[modes]
    ax[1].bar(modes, d.between, color="#c8ccd2", label="between-subject (static)")
    ax[1].bar(modes, d.mean_dynamics, bottom=d.between, color="#f2a541",
              label="mean dynamics μ(t)")
    ax[1].bar(modes, d.subject_dynamics, bottom=d.between + d.mean_dynamics,
              color="#d1495b", label="subject-specific dynamics")
    for i, m in enumerate(modes):
        ax[1].text(i, 1.01, f"{d.within.iloc[i]:.0%} temporal", ha="center", fontsize=9)
    ax[1].set_ylim(0, 1.1)
    ax[1].set_ylabel("share of in-mask pixel variance")
    ax[1].set_title("Variance decomposition (systole)")
    ax[1].legend(frameon=False, fontsize=8, loc="lower right")
    ax[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("preps", nargs="+", help="preprocessing output dirs")
    ap.add_argument("--report", default="application/echonet/report")
    ap.add_argument("--n-train", type=int, default=2000)
    ap.add_argument("--k", type=int, nargs="+", default=[10, 20, 50])
    ap.add_argument("--seed", type=int, default=20261123)
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw the summary figure from existing csv tables")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.report, "figures"), exist_ok=True)
    if a.plot_only:
        ef_figure(pd.read_csv(os.path.join(a.report, "ef_baselines.csv")),
                  pd.read_csv(os.path.join(a.report, "variance_decomposition.csv")),
                  os.path.join(a.report, "figures", "phase0_summary.png"))
        return

    vrows, erows = [], []
    for prep in a.preps:
        mode = os.path.basename(os.path.normpath(prep)).replace("prep_", "")
        v = pd.read_csv(os.path.join(prep, "videos.csv"))
        v = v[v.exclusion.isna() | (v.exclusion == "")]
        mask = np.load(os.path.join(prep, "mask.npy"))
        tr = v[v.Split == "TRAIN"].sample(min(a.n_train, (v.Split == "TRAIN").sum()),
                                          random_state=a.seed)
        va, te = v[v.Split == "VAL"], v[v.Split == "TEST"]
        print(f"[{mode}] train={len(tr)} val={len(va)} test={len(te)} "
              f"mask={mask.sum()}", flush=True)
        X = {s: load_split(prep, d.id.values, mask)
             for s, d in (("tr", tr), ("va", va), ("te", te))}

        for norm in (False, True):
            Xn = X["tr"]
            if norm:
                Xn = Xn / Xn.mean(axis=(1, 2), keepdims=True)
            vrows.append(dict(alignment=mode, normalised=norm,
                              mask_pixels=int(mask.sum()), **var_decomp(Xn)))
        maps_figure(X["tr"], mask,
                    f"alignment = {mode}  (n = {len(tr)} TRAIN videos, systole)",
                    os.path.join(a.report, "figures", f"variance_maps_{mode}.png"))

        F = {s: features(X[s]) for s in X}
        y = dict(tr=tr.EF.values, va=va.EF.values, te=te.EF.values)
        for name in F["tr"]:
            for k in a.k:
                r = ridge_eval(F["tr"][name], y["tr"], F["va"][name], y["va"],
                               F["te"][name], y["te"], k)
                erows.append(dict(alignment=mode, features=name, k=k, **r))
                print(f"  {name:22s} k={k:3d} R2={r['R2']:.3f} "
                      f"MAE={r['MAE']:.2f} AUC={r['AUC_EF40']:.3f}", flush=True)

    vt = pd.DataFrame(vrows)
    et = pd.DataFrame(erows)
    vt.to_csv(os.path.join(a.report, "variance_decomposition.csv"), index=False)
    et.to_csv(os.path.join(a.report, "ef_baselines.csv"), index=False)
    with open(os.path.join(a.report, "tables.md"), "w") as fh:
        fh.write("## Variance decomposition (in-mask, registered systole)\n\n")
        fh.write(vt.round(3).to_markdown(index=False) + "\n\n")
        fh.write("## EF prediction from PCA baselines (TEST)\n\n")
        fh.write(et.round(3).to_markdown(index=False) + "\n")
    ef_figure(et, vt, os.path.join(a.report, "figures", "phase0_summary.png"))


if __name__ == "__main__":
    main()
