#pragma once

#include <eigen3/Eigen/Dense>

namespace spinchain{

using namespace Eigen;


/*! \brief A normalized matrix split into left * right, keeping only the largest Schmidt values
 */
struct SchmidtSplit {
    MatrixXcd left;
    MatrixXcd right;
    VectorXd schmidtValues; // kept Schmidt values in descending order, with norm 1
    double discardedWeight; // fraction of the squared norm that was dropped
};


/*! \brief Truncated Schmidt decomposition of a two-site wavefunction
 *
 *  If matrix = U S V^dagger, this returns left = U and right = S V^dagger when
 *  centerRight is true, or left = U S and right = V^dagger otherwise, restricted to
 *  the largest Schmidt values. The factor holding S is normalized.
 *
 *  Schmidt values are dropped from the smallest up, as long as at most maxKeep
 *  remain or the total discarded weight stays within cutoff. Values at the level
 *  of rounding noise are always dropped, and at least one is always kept.
 *
 *  This diagonalizes the reduced density matrix rather than calling an SVD routine.
 *  Eigen's BDCSVD reads out of bounds in perturbCol0 (Eigen 3.4.0) on the nearly rank
 *  deficient matrices that DMRG produces, and JacobiSVD is far too slow at large bond
 *  dimension. The cost is that Schmidt values below about 1e-8 are not resolved, but
 *  their weight is below 1e-16. The factor holding S is computed by projecting matrix
 *  onto the kept vectors, so no Schmidt value is ever divided by.
 *
 *  Throws std::runtime_error if the matrix is zero or not finite.
 */
SchmidtSplit schmidtSplit(const MatrixXcd& matrix, Index maxKeep, double cutoff, bool centerRight);

} // namespace spinchain
