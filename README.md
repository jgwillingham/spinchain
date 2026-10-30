# SpinChain: Quantum Spin Chains in C++

This project serves two purposes:

First, it provides a basic library for studying quantum spin chains using exact diagonalization.

Second, it serves as an example for any of my future C++ projects. It is meant to illustrate good styling and design practices.


## Arbitrary Spin Heisenberg Model

The code implements an arbitrary spin Heisenberg XYZ model on a chain where the couplings $J_x$, $J_y$, and $J_z$ can be set however you like. The Hamiltonian is

$$ H = \sum_{i=0}^{L-1} J_x s_i^xs_{i+1}^x + J_y s_i^y s_{i+1}^y + J_zs_i^z s_{i+1}^z  - \sum_{i=0}^{L-1} \boldsymbol{h}\cdot \boldsymbol{s}_i$$

where the spin operators are constructed for arbitrary spin using the standard Lie algebra relations. The external field $\boldsymbol{h}=(h_x,h_y,h_z)$ defaults to zero.

Note that this library builds the exact Hamiltonian, so the state space dimension grows exponentially with chain length. Any practical calculation with this must keep chains very short.
