#!/usr/bin/env python3
"""Frozen RETFound-OCT encoder: one 196 x 1024 patch matrix per eye-visit."""
import argparse
import csv
from functools import partial
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
DEPS = ROOT/'application/olives/output/python-deps'
if DEPS.exists():
    sys.path.insert(0, str(DEPS))
import numpy as np
from PIL import Image
import torch
from torch import nn
import timm
from timm.models.vision_transformer import VisionTransformer
import torchvision
from torchvision import transforms

CHECKPOINT_URL = 'https://drive.usercontent.google.com/download?id=1m6s7QYkjyjJDlpEuXm7Xp3PmjN-elfW2&export=download&confirm=t'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def evaluation_transform():
    # Exact deterministic evaluation recipe in openmedlab/RETFound_MAE/util/datasets.py.
    return transforms.Compose([
        transforms.Resize(256, interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.CenterCrop(224), transforms.ToTensor(),
        transforms.Normalize((.485, .456, .406), (.229, .224, .225)),
    ])


def load_encoder(checkpoint):
    # Same ViT architecture and parameter names as the official model. Modern timm
    # supplies the implementation, avoiding the original Python 3.7/timm 0.3 stack.
    model = VisionTransformer(img_size=224, patch_size=16, in_chans=3, num_classes=0,
                              embed_dim=1024, depth=24, num_heads=16, mlp_ratio=4,
                              qkv_bias=True, norm_layer=partial(nn.LayerNorm, eps=1e-6))
    with torch.serialization.safe_globals([argparse.Namespace]):
        checkpoint_data = torch.load(checkpoint, map_location='cpu', weights_only=True, mmap=True)
    state = checkpoint_data['model'] if 'model' in checkpoint_data else checkpoint_data
    expected = model.state_dict()
    missing = sorted(set(expected)-set(state))
    mismatched = [k for k in expected if k in state and expected[k].shape != state[k].shape]
    ignored = sorted(set(state)-set(expected))
    # Pretraining decoder and optimizer are unused. Never tolerate random encoder weights.
    unexpected_encoder = [k for k in ignored if not k.startswith(('decoder_', 'mask_token', 'head.'))]
    if missing or mismatched or unexpected_encoder:
        raise ValueError(dict(missing=missing, mismatched=mismatched, unexpected_encoder=unexpected_encoder))
    model.load_state_dict({k:state[k] for k in expected}, strict=True)
    model.eval().requires_grad_(False)
    return model, dict(encoder_keys_loaded=len(expected), ignored_pretraining_keys=ignored)


@torch.inference_mode()
def encode(model, batch):
    # No random MAE masking, no decoder, no classifier. CLS participates in attention
    # but is removed from the saved spatial matrix after final encoder normalization.
    tokens = model.forward_features(batch)
    if tokens.shape[1:] != (197, 1024):
        raise ValueError(f'Unexpected token shape {tokens.shape}')
    return tokens[:, 1:, :], tokens[:, 0, :]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir', type=Path, default=ROOT/'application/olives/output/retfound')
    ap.add_argument('--checkpoint', type=Path, default=ROOT/'application/olives/output/checkpoints/RETFound_oct_weights.pth')
    ap.add_argument('--batch-size', type=int, default=8)
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--limit', type=int, help='Smoke test only; saved under smoke_N, never overwrites full output')
    args = ap.parse_args()
    if args.batch_size < 1 or (args.limit is not None and args.limit < 1):
        ap.error('batch-size and limit must be positive')
    torch.set_num_threads(min(4, int(os.environ.get('SLURM_CPUS_PER_TASK', 4))))
    torch.manual_seed(20260905)
    if args.device.startswith('cuda') and not torch.cuda.is_available():
        raise RuntimeError('CUDA unavailable: run on an allocated GPU node, or explicitly use --device cpu')
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    manifest = args.input_dir/'manifest.csv'
    with manifest.open() as f:
        rows = list(csv.DictReader(f))
    if [int(r['feature_row']) for r in rows] != list(range(len(rows))):
        raise ValueError('Manifest feature_row must be contiguous and zero-based')
    if args.limit:
        rows = rows[:args.limit]
    out = args.input_dir/(f'smoke_{args.limit}' if args.limit else 'features')
    out.mkdir(parents=True, exist_ok=True)
    completion = out/'extraction.json'
    if completion.exists():
        raise FileExistsError(f'Completed features already exist: {out}; use a new input directory for another run')
    print('Loading pretrained encoder', flush=True)
    model, loading = load_encoder(args.checkpoint)
    model.to(args.device)
    transform = evaluation_transform()
    n = len(rows)
    paths = {k:out/f'{k}.npy.partial' for k in ['patch_tokens', 'pooled', 'cls']}
    arrays = {
        'patch_tokens':np.lib.format.open_memmap(paths['patch_tokens'], mode='w+', dtype='float32', shape=(n,196,1024)),
        'pooled':np.lib.format.open_memmap(paths['pooled'], mode='w+', dtype='float32', shape=(n,1024)),
        'cls':np.lib.format.open_memmap(paths['cls'], mode='w+', dtype='float32', shape=(n,1024)),
    }
    start = time.time()
    repeat_max = None
    for first in range(0,n,args.batch_size):
        imgs = []
        for r in rows[first:first+args.batch_size]:
            with Image.open(ROOT/'data/OLIVES'/r['relative_path']) as im:
                grayscale = im.convert('L')
                if hashlib.sha256(grayscale.tobytes()).hexdigest() != r['pixel_sha256']:
                    raise ValueError(f'Source image changed: {r["relative_path"]}')
                imgs.append(transform(im.convert('RGB')))
        batch = torch.stack(imgs).to(args.device)
        patch, cls = encode(model, batch)
        if first == 0:
            again, _ = encode(model, batch)
            repeat_max = float((again-patch).abs().max())
            if repeat_max > 1e-5:
                raise ValueError(f'Nonrepeatable features: {repeat_max}')
            # Verify the library forward against the official unmasked encoder sequence.
            z = model.patch_embed(batch)
            z = torch.cat((model.cls_token.expand(len(batch),-1,-1), z), dim=1)
            z = model.pos_drop(z + model.pos_embed)
            for block in model.blocks:
                z = block(z)
            z = model.norm(z)
            reference_max = float((z[:,1:]-patch).abs().max())
            if reference_max > 1e-5:
                raise ValueError(f'Encoder differs from official sequence: {reference_max}')
        values = {'patch_tokens':patch, 'pooled':patch.mean(dim=1), 'cls':cls}
        for key,value in values.items():
            a = value.cpu().numpy()
            if not np.isfinite(a).all():
                raise ValueError(f'Nonfinite {key} at row {first}')
            arrays[key][first:first+len(imgs)] = a
        print(f'{min(first+args.batch_size,n)}/{n} images, {time.time()-start:.1f}s', flush=True)
    for a in arrays.values():
        a.flush()
    del arrays
    for key,path in paths.items():
        path.replace(out/f'{key}.npy')
    # A raw binary companion lets the R MatFPCA implementation consume the
    # matrices without depending on reticulate or a particular NumPy version.
    # The order is eye-visit, patch, feature, with feature varying fastest.
    patch_map = np.load(out/'patch_tokens.npy', mmap_mode='r')
    patch_map.transpose(0, 1, 2).astype('<f4', copy=False).tofile(out/'patch_tokens.bin')
    record = dict(model='RETFound OCT MAE ViT-L/16', frozen=True,
                  checkpoint_url=CHECKPOINT_URL, checkpoint_sha256=sha256(args.checkpoint),
                  source_repository='https://github.com/openmedlab/RETFound_MAE',
                  manifest_sha256=sha256(manifest), n_images=n, patch_shape=[n,196,1024], dtype='float32',
                  preprocessing='PIL RGB; bicubic shorter-side resize 256; center crop 224; ImageNet mean/std',
                  laterality_flip=False, axial_alignment=False, random_augmentation=False,
                  tokens='final encoder LayerNorm; patch_tokens excludes CLS; cls.npy stores the global CLS token; patch order is raster (14 rows x 14 columns)',
                  pooled='mean of final normalized patch tokens; no task-specific fc_norm or classifier',
                  matrix_axes=['eye_visit','spatial_patch','learned_feature'],
                  repeat_max_abs_difference=repeat_max, official_sequence_max_abs_difference=reference_max,
                  encoder_load=loading, torch=torch.__version__, torchvision=torchvision.__version__,
                  timm=timm.__version__, python=platform.python_version(), device=args.device,
                  gpu=torch.cuda.get_device_name() if args.device.startswith('cuda') else None,
                  seconds=time.time()-start)
    completion.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)


if __name__ == '__main__':
    main()
