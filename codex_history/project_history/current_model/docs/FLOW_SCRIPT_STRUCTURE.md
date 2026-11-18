# Flow script structure

Every flow Lua script should use the following top-level blocks, in order:

1. `Model summary`
2. `Dependencies`
3. `Parameters`
4. `Discretization and boundary conditions`
5. `Mesh and approximation space`
6. `Solver`
7. `Time discretization and initial values`
8. `Output`

Each block begins with the same dashed separator and contains a short comment
explaining its purpose. Boundary and coefficient functions are declared at the
top of the discretization block, before they are registered with UG4.

## Runtime input contract

Flow scripts expose only these command-line parameters:

- `-numRefs`: number of mesh refinements
- `-output`: VTK output prefix

Mesh paths, physical coefficients, boundary values, discretization options,
and solver settings belong in the script's `Parameters` block. This keeps a
simulation configuration reproducible without reconstructing a terminal
command. An MPI launch remains supported, for example:

```sh
mpirun -np 8 ugshell -ex code/flow/example.lua -numRefs 2 -output Results/run
```

Steady scripts still include the `Time discretization and initial values`
block; its summary states that no time integrator is used and identifies how
the initial algebraic state is constructed.
