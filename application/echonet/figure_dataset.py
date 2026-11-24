#!/usr/bin/env python
"""Dataset explainer figure (PLAN.md Section 3, item 1).

Top: ED and ES frames for three patients (low, normal, high EF) with the
expert LV tracing overlaid. Bottom: mean grey level inside the ED LV outline
over the whole clip (the blood pool is dark, so the curve rises as the LV
empties), with ED/ES marked, and the systolic window mapped to phase [0, 1].
"""
import argparse
import os
import sys

import av
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.path import Path

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from prepare_echonet import read_landmarks  # noqa: E402


def outline(vt, vid, frame):
    g = vt[(vt.id == vid) & (vt.Frame == frame)].iloc[1:]
    return np.r_[g[["X1", "Y1"]].values, g[["X2", "Y2"]].values[::-1]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/EchoNet-Dynamic")
    ap.add_argument("--out", default="application/echonet/report/figures")
    ap.add_argument("--seed", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fl = pd.read_csv(os.path.join(a.data, "FileList.csv")).rename(columns={"FileName": "id"})
    vt = pd.read_csv(os.path.join(a.data, "VolumeTracings.csv"))
    vt["id"] = vt.FileName.str[:-4]
    lm = read_landmarks(os.path.join(a.data, "VolumeTracings.csv"))
    v = fl.merge(lm, on="id")
    v = v[(v.ES > v.ED) & (v.FrameHeight == 112) & (v.FPS == 50)
          & ((v.ES - v.ED) / v.FPS).between(0.2, 0.45)]
    rng = np.random.default_rng(a.seed)
    picks = []
    for lo, hi, lab in [(15, 30, "low EF"), (55, 62, "normal EF"), (72, 85, "high EF")]:
        c = v[v.EF.between(lo, hi)]
        picks.append((c.iloc[rng.integers(len(c))], lab))

    fig = plt.figure(figsize=(13, 7))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 0.85])
    for j, (r, lab) in enumerate(picks):
        with av.open(os.path.join(a.data, "Videos", r.id + ".avi")) as c:
            g = np.stack([f.to_ndarray(format="rgb24") for f in c.decode(video=0)]).mean(-1)
        ED, ES = int(r.ED), int(r.ES)
        for k, (fr, nm, col) in enumerate([(ED, "ED", "tab:blue"), (ES, "ES", "tab:red")]):
            ax = fig.add_subplot(gs[0, 2 * j + k])
            ax.imshow(g[fr], cmap="gray", vmin=0, vmax=255)
            P = outline(vt, r.id, fr)
            ax.plot(*np.r_[P, P[:1]].T, color=col, lw=1.2)
            ax.set_title(f"{lab} ({r.EF:.0f}%)\n{nm}, frame {fr}" if k == 0
                         else f"\n{nm}, frame {fr}", fontsize=9)
            ax.axis("off")
        # LV intensity curve inside ED outline
        P = outline(vt, r.id, ED)
        yy, xx = np.mgrid[:112, :112]
        inside = Path(P).contains_points(np.c_[xx.ravel(), yy.ravel()]).reshape(112, 112)
        curve = g[:, inside].mean(1)
        t = np.arange(len(g)) / r.FPS
        ax = fig.add_subplot(gs[1, 2 * j:2 * j + 2])
        ax.plot(t, curve, color="0.3", lw=1)
        ax.axvspan(ED / r.FPS, ES / r.FPS, color="tab:orange", alpha=0.2,
                   label="systole → phase [0, 1]")
        ax.axvline(ED / r.FPS, color="tab:blue", lw=1.2, label="ED (phase 0)")
        ax.axvline(ES / r.FPS, color="tab:red", lw=1.2, label="ES (phase 1)")
        ax.set_xlabel("clock time (s)")
        if j == 0:
            ax.set_ylabel("mean grey inside ED LV outline")
            ax.legend(fontsize=8, loc="upper right")
        ax.set_title(f"{r.FPS:.0f} fps, {len(g)} frames; {ES - ED + 1} systolic frames",
                     fontsize=9)
    fig.suptitle("EchoNet-Dynamic: apical four-chamber echo videos with expert ED/ES landmarks")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "dataset_explainer.png"), dpi=90)
    print("wrote", os.path.join(a.out, "dataset_explainer.png"))


if __name__ == "__main__":
    main()
