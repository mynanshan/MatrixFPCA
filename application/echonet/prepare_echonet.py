#!/usr/bin/env python
"""EchoNet-Dynamic preprocessing for MatFPCA (PLAN.md Sections 2.1, 2.2).

For every video this script
  1. reads the expert ED/ES landmarks from VolumeTracings.csv,
  2. applies the eligibility rules and records an exclusion reason
     (frame size, tracings, ES before ED, ED->ES gap, implausible long axis,
     colour-map display),
  3. decodes the clip with PyAV and converts it to grey (mean of RGB),
  4. removes overlays outside the ultrasound fan (depth ticks, apex marker),
  5. optionally aligns the left ventricle with one transform per video taken
     from the ED long axis (--align none | rigid | similarity),
  6. keeps the systolic frames ED..ES and assigns cardiac phase
     t = (frame - ED) / (ES - ED), snapped to a 51-point grid,
  7. writes the registered frames as uint8 to <out>/frames/<id>.npz.

A second pass builds the consensus sector mask from TRAIN videos and the
per-video intensity-normalisation scale, then writes
  videos.csv      one row per video (landmarks, QC, exclusion reason, scale)
  manifest.csv    one row per retained frame
  qc_summary.csv  exclusion counts
  mask.npy/.png   consensus analysis mask
  config.json     run settings and the canonical LV landmarks
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import av
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import Delaunay

GRID_STEP = 0.02
N_GRID = 51
SIZE = 112


# --------------------------------------------------------------------------
# landmarks
# --------------------------------------------------------------------------
def read_landmarks(tracings_csv):
    """ED/ES frames and the ED long axis (apex, base) for each video.

    Each traced frame has the long axis in its first row and perpendicular
    chords after it. Area ~ sum(chord lengths) * axis length / n_chords.
    ED is the traced frame with the larger area.
    """
    vt = pd.read_csv(tracings_csv)
    vt["id"] = vt.FileName.str.replace(r"\.avi$", "", regex=True)
    rows = []
    for (vid, frame), g in vt.groupby(["id", "Frame"], sort=False):
        ax = g.iloc[0]
        ch = g.iloc[1:]
        L = np.hypot(ax.X2 - ax.X1, ax.Y2 - ax.Y1)
        chords = np.hypot(ch.X2 - ch.X1, ch.Y2 - ch.Y1)
        area = chords.sum() * L / max(len(ch), 1)
        # apex is the end nearer the top of the image (smaller y)
        if ax.Y1 <= ax.Y2:
            apex, base = (ax.X1, ax.Y1), (ax.X2, ax.Y2)
        else:
            apex, base = (ax.X2, ax.Y2), (ax.X1, ax.Y1)
        rows.append(dict(id=vid, frame=int(frame), area=area,
                         apex_x=apex[0], apex_y=apex[1],
                         base_x=base[0], base_y=base[1]))
    a = pd.DataFrame(rows)
    out = []
    for vid, g in a.groupby("id", sort=False):
        if len(g) != 2:
            out.append(dict(id=vid, n_traced=len(g)))
            continue
        g = g.sort_values("area")
        es, ed = g.iloc[0], g.iloc[1]
        out.append(dict(id=vid, n_traced=2, ED=int(ed.frame), ES=int(es.frame),
                        area_ED=ed.area, area_ES=es.area,
                        apex_x=ed.apex_x, apex_y=ed.apex_y,
                        base_x=ed.base_x, base_y=ed.base_y))
    return pd.DataFrame(out)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
def lv_transform(apex, base, canon_apex, canon_base, mode):
    """Return (matrix, offset) for ndimage.affine_transform (output -> input),
    in (row, col) coordinates, mapping this video's LV onto the canonical one.

    rigid:       rotate the long axis to the canonical angle and move the LV
                 centre onto the canonical centre (no scaling, keeps size).
    similarity:  additionally scale so the long axis has canonical length.
    """
    z_a, z_b = complex(*apex), complex(*base)          # x + i y
    c_a, c_b = complex(*canon_apex), complex(*canon_base)
    rot = (c_b - c_a) / (z_b - z_a)                     # s * exp(i theta)
    if mode == "rigid":
        rot = rot / abs(rot)
    z_c, c_c = (z_a + z_b) / 2, (c_a + c_b) / 2
    # forward: z_out = rot * (z_in - z_c) + c_c
    # inverse: z_in  = (z_out - c_c) / rot + z_c  = A z_out + b
    A = 1 / rot
    b = z_c - c_c * A
    # in (row=y, col=x): [y_in, x_in] = M [y_out, x_out] + off
    ar, ai = A.real, A.imag
    M = np.array([[ar, ai], [-ai, ar]])
    off = np.array([b.imag, b.real])
    return M, off


def fan_mask(grey_max, thr=8.0):
    """Ultrasound fan of one clip: largest connected bright region after a
    small opening, holes filled. Removes the depth-scale ticks and the apex
    marker, which are thin and detached from the fan."""
    m = grey_max > thr
    m = ndimage.binary_opening(m, structure=np.ones((3, 3)))
    lab, n = ndimage.label(m)
    if n == 0:
        return m
    sizes = ndimage.sum(m, lab, index=np.arange(1, n + 1))
    m = lab == (np.argmax(sizes) + 1)
    # the sector is convex; take the hull so dark (attenuated) regions at the
    # bottom of the fan are not cut out
    pts = np.argwhere(m)
    hull = Delaunay(pts)
    rr, cc = np.mgrid[:m.shape[0], :m.shape[1]]
    inside = hull.find_simplex(np.c_[rr.ravel(), cc.ravel()]) >= 0
    return inside.reshape(m.shape)


# --------------------------------------------------------------------------
# per-video worker
# --------------------------------------------------------------------------
def process_video(args):
    r, video_dir, out_dir, mode, canon = args
    vid = r["id"]
    res = dict(id=vid)
    try:
        with av.open(os.path.join(video_dir, vid + ".avi")) as c:
            rgb = np.stack([f.to_ndarray(format="rgb24")
                            for f in c.decode(video=0)])
    except Exception as e:  # noqa: BLE001
        res["exclusion"] = f"decode_error: {e}"
        return res, None, None
    n_dec = len(rgb)
    res["n_decoded"] = n_dec
    ED, ES = int(r["ED"]), int(r["ES"])
    if ES >= n_dec:
        res["exclusion"] = "landmark_beyond_clip"
        return res, None, None

    rgb = rgb.astype(np.int16)
    chroma = rgb.max(-1) - rgb.min(-1)
    grey = rgb.mean(-1).astype(np.float32)               # 0..255
    fan = fan_mask(grey.max(0))
    res["fan_pixels"] = int(fan.sum())
    res["colour_frac"] = float((chroma[:, fan] > 20).mean())  # in-fan colour
    res["overlay_frac"] = float((chroma[:, ~fan] > 20).mean())

    sys_frames = grey[ED:ES + 1] * fan                    # overlays removed
    if mode != "none":
        M, off = lv_transform((r["apex_x"], r["apex_y"]),
                              (r["base_x"], r["base_y"]),
                              canon["apex"], canon["base"], mode)
        sys_frames = np.stack([ndimage.affine_transform(
            f, M, offset=off, order=1, mode="constant", cval=0.0)
            for f in sys_frames])
        fan_w = ndimage.affine_transform(fan.astype(np.float32), M,
                                         offset=off, order=1) > 0.99
        sys_frames *= fan_w
    else:
        fan_w = fan

    # cardiac phase and grid snapping
    idx = np.arange(ED, ES + 1)
    phase = (idx - ED) / (ES - ED)
    gidx = np.rint(phase / GRID_STEP).astype(int)
    dist = np.abs(phase - gidx * GRID_STEP)
    keep = np.zeros(len(idx), bool)
    for g in np.unique(gidx):                              # closest per grid pt
        w = np.where(gidx == g)[0]
        keep[w[np.argmin(dist[w])]] = True
    frames = np.clip(np.rint(sys_frames[keep]), 0, 255).astype(np.uint8)

    np.savez_compressed(os.path.join(out_dir, "frames", vid + ".npz"),
                        frames=frames, frame_idx=idx[keep],
                        phase=phase[keep], grid_idx=gidx[keep],
                        fan=fan_w)
    res.update(n_retained=int(keep.sum()), n_systolic=len(idx),
               exclusion="")
    frame_rows = pd.DataFrame(dict(id=vid, frame_idx=idx, phase=phase,
                                   grid_idx=gidx, retained=keep))
    return res, frame_rows, fan_w


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--data", default="data/EchoNet-Dynamic")
    ap.add_argument("--out", default="application/echonet/output/prep_none")
    ap.add_argument("--align", choices=["none", "rigid", "similarity"],
                    default="none")
    ap.add_argument("--gap-min", type=float, default=0.10, help="seconds")
    ap.add_argument("--gap-max", type=float, default=0.80, help="seconds")
    ap.add_argument("--mask-coverage", type=float, default=0.90,
                    help="pixel kept if inside the (aligned) fan of at least "
                         "this fraction of TRAIN videos")
    ap.add_argument("--max-colour", type=float, default=0.05,
                    help="exclude videos whose in-fan colour fraction exceeds "
                         "this (colour-map display)")
    ap.add_argument("--limit", type=int, default=0,
                    help="process only a random subset of eligible videos "
                         "(per split, for pilots); 0 = all")
    ap.add_argument("--seed", type=int, default=20261123)
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    a = ap.parse_args()

    t0 = time.time()
    os.makedirs(os.path.join(a.out, "frames"), exist_ok=True)
    fl = pd.read_csv(os.path.join(a.data, "FileList.csv"))
    fl = fl.rename(columns={"FileName": "id"})
    fl["id"] = fl.id.str.replace(r"\.avi$", "", regex=True)
    lm = read_landmarks(os.path.join(a.data, "VolumeTracings.csv"))
    v = fl.merge(lm, on="id", how="left")
    v["gap_s"] = (v.ES - v.ED) / v.FPS

    # eligibility, first matching reason wins
    reason = pd.Series("", index=v.index)
    def flag(cond, why):
        sel = (reason == "") & cond
        reason[sel] = why
    flag((v.FrameHeight != SIZE) | (v.FrameWidth != SIZE), "frame_size")
    flag(v.n_traced.isna(), "no_tracing")
    flag(v.n_traced.fillna(2) != 2, "tracing_count")
    flag(v.ES <= v.ED, "es_before_ed")
    flag(v.gap_s < a.gap_min, "gap_short")
    flag(v.gap_s > a.gap_max, "gap_long")
    # implausible ED long axis (bad tracing): robust z > 4 on length or tilt,
    # judged against the TRAIN videos that passed the rules above
    v["axis_len"] = np.hypot(v.base_x - v.apex_x, v.base_y - v.apex_y)
    v["axis_tilt"] = np.degrees(np.arctan2(v.base_x - v.apex_x,
                                           v.base_y - v.apex_y))
    ref = v[(reason == "") & (v.Split == "TRAIN")]
    for col in ("axis_len", "axis_tilt"):
        med = ref[col].median()
        mad = 1.4826 * (ref[col] - med).abs().median()
        flag((v[col] - med).abs() > 4 * mad, "landmark_outlier")
    v["exclusion"] = reason
    elig = v[reason == ""]
    if a.limit:
        elig = pd.concat([g.sample(min(a.limit, len(g)), random_state=a.seed)
                          for _, g in elig.groupby("Split")])

    # canonical ED long axis = median over eligible TRAIN videos
    tr = v[(reason == "") & (v.Split == "TRAIN")]
    canon = dict(apex=(float(tr.apex_x.median()), float(tr.apex_y.median())),
                 base=(float(tr.base_x.median()), float(tr.base_y.median())))

    video_dir = os.path.join(a.data, "Videos")
    jobs = [(r, video_dir, a.out, a.align, canon)
            for r in elig.to_dict("records")]
    print(f"[prep] align={a.align} videos={len(jobs)} workers={a.workers}",
          flush=True)
    results, frame_tabs, fans = [], [], {}
    with Pool(a.workers) as pool:
        for k, (res, fr, fan_w) in enumerate(
                pool.imap_unordered(process_video, jobs, chunksize=8)):
            results.append(res)
            if fr is not None:
                frame_tabs.append(fr)
                fans[res["id"]] = fan_w
            if (k + 1) % 1000 == 0:
                print(f"[prep] {k + 1}/{len(jobs)} {time.time() - t0:.0f}s",
                      flush=True)

    res = pd.DataFrame(results).set_index("id")
    v = v.set_index("id")
    for c in res.columns:
        if c == "exclusion":
            sel = res.index
            v.loc[sel, "exclusion"] = res["exclusion"].fillna("")
        else:
            v.loc[res.index, c] = res[c]
    if a.limit:
        v.loc[(v.exclusion == "") & ~v.index.isin(res.index),
              "exclusion"] = "not_sampled"

    # videos displayed with a colour map (blue tint) instead of grey scale:
    # the RGB mean is a non-linear transform of the echo intensity there
    v.loc[(v.exclusion == "") & (v.colour_frac > a.max_colour),
          "exclusion"] = "colour_tint"

    # consensus mask from retained TRAIN videos
    tr_ids = v.index[(v.exclusion == "") & (v.Split == "TRAIN")]
    coverage = np.mean([fans[i] for i in tr_ids], axis=0)
    mask = coverage >= a.mask_coverage
    np.save(os.path.join(a.out, "mask.npy"), mask)
    np.save(os.path.join(a.out, "mask_coverage.npy"), coverage)
    try:
        from PIL import Image
        Image.fromarray((mask * 255).astype(np.uint8)).save(
            os.path.join(a.out, "mask.png"))
    except Exception:  # noqa: BLE001
        pass

    # intensity normalisation: mean in-mask grey over retained frames (0..1)
    kept = v.index[v.exclusion == ""]
    scale = {}
    for vid in kept:
        f = np.load(os.path.join(a.out, "frames", vid + ".npz"))["frames"]
        scale[vid] = float(f[:, mask].mean() / 255.0)
    v["norm_scale"] = pd.Series(scale)
    v.loc[(v.exclusion == "") & (v.norm_scale < 0.02), "exclusion"] = "blank_in_mask"
    kept = v.index[v.exclusion == ""]
    v = v.reset_index()
    v.to_csv(os.path.join(a.out, "videos.csv"), index=False)

    man = pd.concat(frame_tabs, ignore_index=True)
    man = man.merge(v[["id", "Split", "EF", "ESV", "EDV", "FPS", "ED", "ES",
                       "colour_frac", "norm_scale", "exclusion"]], on="id")
    man = man[man.exclusion == ""].drop(columns="exclusion")
    man.to_csv(os.path.join(a.out, "manifest.csv"), index=False)

    qc = (v.exclusion.replace("", "retained").value_counts()
          .rename_axis("status").reset_index(name="n_videos"))
    by_split = pd.crosstab(v.exclusion.replace("", "retained"), v.Split)
    qc = qc.merge(by_split, left_on="status", right_index=True)
    qc.to_csv(os.path.join(a.out, "qc_summary.csv"), index=False)

    cfg = dict(vars(a), canonical=canon, mask_pixels=int(mask.sum()),
               n_retained=int(len(kept)), seconds=round(time.time() - t0, 1))
    with open(os.path.join(a.out, "config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2)
    print(qc.to_string(index=False))
    print(f"[prep] mask pixels={mask.sum()} done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
