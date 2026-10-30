#pragma once

#include <MPO.hpp>
#include <MPS.hpp>
#include <Tensors.hpp>

namespace spinchain{


/*  An environment is the Hamiltonian contracted with the state over part of the
 *  chain, with indices (bra, mpo, ket):
 *
 *      .--- bra                    bra ---.
 *      |                                  |
 *      L--- mpo                    mpo ---R
 *      |                                  |
 *      .--- ket                    ket ---.
 */

// the environment of an empty part of the chain
Tensor3 emptyEnvironment();

// extends a left environment across a site with MPS tensor A and MPO tensor W
Tensor3 growLeft(const Tensor3& left, const Tensor3& A, const Tensor4& W);

// extends a right environment across a site with MPS tensor A and MPO tensor W
Tensor3 growRight(const Tensor3& right, const Tensor3& A, const Tensor4& W);

/*! \brief Applies the effective two-site Hamiltonian to theta without ever forming it
 *
 *      .---     ---.
 *      |    W1 W2  |
 *      L----|--|---R
 *      |    |  |   |
 *      .--- theta -.
 *
 *  theta has indices (bondLeft, physical1, physical2, bondRight), as does the result.
 */
Tensor4 applyTwoSiteHamiltonian(const Tensor3& left, const Tensor4& W1, const Tensor4& W2,
                                const Tensor3& right, const Tensor4& theta);

// returns <psi|H|psi> for a normalized state
double energy(const MPS& state, const MPO& hamiltonian);

} // namespace spinchain
