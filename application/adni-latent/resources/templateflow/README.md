# TemplateFlow atlas used by the BrainIAC concept audit

`tpl-MNI152NLin2009cAsym_res-01_seg-aseg_dseg.nii.gz` is the official
TemplateFlow MNI152NLin2009cAsym 1 mm `aseg` segmentation. It was downloaded
from:

<https://templateflow.s3.amazonaws.com/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_seg-aseg_dseg.nii.gz>

SHA-256:

```text
0bb033414fbbd6df626d063b73801a14c71fb33b83204585f2aa1fbb4aef60dd
```

The concept audit records the exact integer labels used for each bilateral
region in its `run_manifest.json`. The atlas defines regions only; all
visit-specific physical measurements come from the independently registered
ADNI log-Jacobian maps.
