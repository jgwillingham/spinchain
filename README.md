# SpinChain: Quantum Spin Chains in C++

This project serves two purposes:

First, it provides a basic library for studying quantum spin chains using exact diagonalization and density matrix renormalization group (DMRG).

Second, it serves as an example for any of my future C++ projects. It is meant to illustrate good styling and design practices.


## Arbitrary Spin Heisenberg Model

The code implements an arbitrary spin Heisenberg XYZ model on a chain where the couplings $J_x$, $J_y$, and $J_z$ can be set however you like. The Hamiltonian is

$$ H = \sum_{i=0}^{L-1} J_x s_i^xs_{i+1}^x + J_y s_i^y s_{i+1}^y + J_zs_i^z s_{i+1}^z  - h_x s_i^x - h_y s_i^y - h_z s_i^z$$

where the spin operators are constructed for arbitrary spin using the standard Lie algebra relations. The external field $\boldsymbol{h}=(h_x,h_y,h_z)$ defaults to zero.


## Solvers

| Solver | Finds | Chain length | Accuracy |
| --- | --- | --- | --- |
| `ExactDiagonalization` | whole spectrum | very short: the state space grows exponentially | exact |
| `DMRG` | ground state | long | controlled by the bond dimension $\chi$ |

A simple example program which uses DMRG to find the ground state magnetization of an isotropic spin-1/2 Heisenberg chain of length 100 is:
```cpp
SpinChain model(100, {1.0, 1.0, 1.0}, 0.5);

DMRGOptions options;
options.chiMax = 64;

DMRG solver(model, options);
SolveResult result = solver.solve();
Vector3cd magnetization = solver.magnetization();
```

To add another method, derive from `Solver` and implement `solve` and `localExpectations`.

### DMRG

`DMRG` is a two-site finite DMRG for open chains. It works for any couplings and field, since it uses dense tensors and does not rely on any symmetry of the model.

Accuracy is set by `chiMax`. A chain with a gap converges at a small $\chi$, while a critical chain such as the $S=1/2$ Heisenberg chain needs more as it gets longer. A run with too small a $\chi$ can report that it converged and still be above the true ground energy, so compare a few values of $\chi$ and look at the truncation error. The same is true of the two-site Lanczos solves, where `lanczosIterations` being too small slows the convergence of observables much more than the energy.

If the ground state is degenerate, as it is for an open $S=1$ chain at long lengths, DMRG returns some state in that manifold and only the energy is meaningful.

## Building

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Or use Docker, where `docker build --target test .` builds and runs the tests.

The tests rely on `assert`, so the default build type is Release without `-DNDEBUG`. DMRG runs about twice as fast with `-DCMAKE_CXX_FLAGS_RELEASE="-O2 -march=native"`, but the binary then only runs on similar CPUs.
