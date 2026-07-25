# MatrixFPCA

## Eigenfunction truncation

`trunc.matfeigen()` retains the full row and column eigensystems by default and
applies its PVE threshold only to the nonnegative temporal eigenvalues:

```r
eig_trunc <- trunc.matfeigen(eig, pve = 0.95)
```

The earlier greedy rule, which can trim the row, column, and temporal
eigensystems together, remains available with `pve_rule = "joint"`:

```r
eig_trunc <- trunc.matfeigen(eig, pve = 0.95, pve_rule = "joint")
```

Supplying all three ranks (`I`, `J`, and `K`) bypasses PVE selection.
