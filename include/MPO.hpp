#pragma once

#include <vector>
#include <SpinChain.hpp>
#include <Tensors.hpp>

namespace spinchain{


/*! \brief Matrix product operator for the Hamiltonian of a SpinChain
 *
 *  The Hamiltonian is written as a product of site tensors W[i], with indices
 *  (mpoLeft, mpoRight, physicalOut, physicalIn). With one row for each of the n
 *  bond terms in SpinChain::bondTerms(), the matrix of operators on a bulk site is
 *
 *      [ I  left_1  ...  left_n  onsite        ]
 *      [ 0  0       ...  0       J_1 right_1   ]
 *      [ .  .            .       .             ]
 *      [ 0  0       ...  0       J_n right_n   ]
 *      [ 0  0       ...  0       I             ]
 *
 *  This is a small state machine read from left to right: the first row means no
 *  term has started, the last row means a term has finished, and the rows in
 *  between mean a bond term was started on the previous site. The first site
 *  keeps only the first row and the last site keeps only the last column.
 */
class MPO {
    public:
        explicit MPO(const SpinChain& model);

        std::vector<Tensor4> W;
};

} // namespace spinchain
