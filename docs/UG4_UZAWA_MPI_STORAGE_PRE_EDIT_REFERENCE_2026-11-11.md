# UG4 Uzawa MPI-storage pre-edit reference — 2026-11-11

This is the rollback reference made immediately before testing consistent MPI
storage for the two Uzawa cross-block matrix products.

- Repository: `/Users/alhilaba/UG4_promesh_ogrid/ugcore`
- `HEAD`: `46b9b57a8c8398939bf1b6a4661396c863f6a723`
- File: `ugbase/lib_disc/operator/linear_operator/uzawa/uzawa.h`
- Pre-edit SHA-256:
  `17a735273e78ab330f17f3b4b5e639fb8123ec01f0b85c8b8a1c99de6423b42d`
- Pre-edit Git status for this file: clean (no diff from `HEAD`)

Therefore the exact pre-edit file can be recovered without affecting unrelated
UG4 work with:

```text
git show 46b9b57a8c8398939bf1b6a4661396c863f6a723:ugbase/lib_disc/operator/linear_operator/uzawa/uzawa.h
```

Do not reset or check out the complete UG4-core repository because it contains
other unrelated local modifications.

## Suspected defect under test

In `UzawaBase::step`, the forward and Schur iterators can return `cRegular` and
`cSchur` in `PST_UNIQUE` storage on MPI runs. Those vectors are subsequently
used as inputs to `B21.apply_sub` and `B12.apply_sub`. Distributed matrix-vector
application requires consistent interface input. The test edit converts each
correction to `PST_CONSISTENT` immediately before its cross-block product and
throws if that conversion fails.
