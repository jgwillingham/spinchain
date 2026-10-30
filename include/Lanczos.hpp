#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <eigen3/Eigen/Dense>

namespace spinchain{

using namespace Eigen;


/*! \brief Outcome of a Lanczos calculation
 */
struct LanczosResult {
    double eigenvalue;      // lowest eigenvalue found
    VectorXcd eigenvector;  // its normalized eigenvector
    size_t iterations;      // number of matrix-vector products used
    bool converged;         // whether the residual norm reached the tolerance
};


/*! \brief Lanczos solver for the lowest eigenpair of a Hermitian operator
 *
 *  The operator is only accessed through applyH, a callable mapping a vector x
 *  to H*x, so H never needs to be stored as a matrix.
 *
 *  Every new Krylov vector is orthogonalized against all previous ones. Without
 *  this, rounding errors destroy their orthogonality and spurious copies of
 *  eigenvalues appear. It costs memory, so keep maxIterations small.
 *
 *  The result stays in the symmetry sector of start: a Krylov space cannot leave
 *  it. If start is an exact eigenvector, then that eigenvector is returned.
 */
template <typename MatVec>
LanczosResult lanczos(const MatVec& applyH, const VectorXcd& start,
                      size_t maxIterations = 30, double tolerance = 1e-10)
{
    double startNorm = start.norm();
    if (startNorm == 0.0 || !std::isfinite(startNorm)){
        throw std::invalid_argument("lanczos requires a nonzero, finite start vector");
    }

    // the Krylov space cannot be larger than the vector space
    size_t maxSteps = std::min(maxIterations, static_cast<size_t>(start.size()));

    std::vector<VectorXcd> basis = {start/startNorm};
    std::vector<double> alpha;  // diagonal of the tridiagonal matrix
    std::vector<double> beta;   // subdiagonal of the tridiagonal matrix

    SelfAdjointEigenSolver<MatrixXd> tridiagonal;
    bool converged = false;

    for (size_t j=0; j<maxSteps; j++){
        VectorXcd w = applyH(basis[j]);
        alpha.push_back(basis[j].dot(w).real());

        // orthogonalize twice: one pass is not enough once the basis is nearly dependent
        for (int pass=0; pass<2; pass++){
            for (const VectorXcd& v : basis){
                w -= v*v.dot(w);
            }
        }
        double nextBeta = w.norm();

        VectorXd diagonal = Map<VectorXd>(alpha.data(), alpha.size());
        VectorXd subdiagonal = Map<VectorXd>(beta.data(), beta.size());
        tridiagonal.computeFromTridiagonal(diagonal, subdiagonal);

        // the residual norm is nextBeta times the last component of the Ritz vector
        double ritzValue = tridiagonal.eigenvalues()(0);
        double residual = nextBeta*std::abs(tridiagonal.eigenvectors()(j,0));

        // a vanishing nextBeta means the Krylov space is invariant: there is nothing left to add
        double breakdown = std::numeric_limits<double>::epsilon()*std::max(1.0, std::abs(ritzValue));
        if (residual <= tolerance || nextBeta <= breakdown){
            converged = true;
            break;
        }

        if (j+1 < maxSteps){
            beta.push_back(nextBeta);
            basis.push_back(w/nextBeta);
        }
    }

    size_t steps = alpha.size();
    VectorXcd eigenvector = VectorXcd::Zero(start.size());
    for (size_t j=0; j<steps; j++){
        eigenvector += tridiagonal.eigenvectors()(j,0)*basis[j];
    }
    eigenvector.normalize();

    return {tridiagonal.eigenvalues()(0), eigenvector, steps, converged};
}

} // namespace spinchain
