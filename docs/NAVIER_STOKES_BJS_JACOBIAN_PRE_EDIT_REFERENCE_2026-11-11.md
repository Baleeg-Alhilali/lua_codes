# Navier--Stokes BJS Jacobian pre-edit reference — 2026-11-11

This is the rollback reference created before testing the missing boundary
shape-function factor in `NavierStokesWSBCFV1::add_jac_A_elem`.

- Repository: `/Users/alhilaba/UG4_promesh_ogrid/plugins/NavierStokes`
- `HEAD`: `8df0ce228312397f10638219ba39383dbf59b4ac`
- File: `incompressible/fv1/bnd/wall_sliding_fv1.cpp`
- Exact pre-edit SHA-256:
  `d71eec2b4b0c8a7f6393eb5c6291608f5392ae1e07649a4c1bb7f4f9fb9447e9`

The file already contained user/project modifications relative to `HEAD`, so
the repository must not be reset wholesale. Immediately before this test, the
tangential BJS Jacobian contribution was:

```cpp
J(d1, bf->node_id(), d2, sh) += parallelVel(d1,d2)
    * m_imSlidingFactor[ip] * bf->volume();
```

To restore that exact pre-test expression, remove only the subsequently added
`* bf->shape(sh)` factor from this contribution. The header and registration
source are not changed by this test.

## Reason for the test

The corresponding defect evaluates the boundary velocity as
`sum_sh u(sh) * bf->shape(sh)`. For the linear BJS term, differentiation with
respect to nodal value `u(sh)` therefore requires the same
`bf->shape(sh)` factor in the Jacobian. Without it, the Jacobian is not the
derivative of the defect and its effective boundary stiffness depends on the
number and type of element shape functions.
