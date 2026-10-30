#include <SchmidtSplit.hpp>
#include <algorithm>
#include <stdexcept>


namespace spinchain{

// density matrix eigenvalues this far below the largest are rounding noise, not physics
static const double NOISE_FLOOR = 1e-15;


SchmidtSplit schmidtSplit(const MatrixXcd& matrix, Index maxKeep, double cutoff, bool centerRight)
{
    // the eigenvectors of the reduced density matrix are the Schmidt vectors of that side
    MatrixXcd densityMatrix;
    if (centerRight){
        densityMatrix = matrix*matrix.adjoint();
    } else {
        densityMatrix = matrix.adjoint()*matrix;
    }

    SelfAdjointEigenSolver<MatrixXcd> solution(densityMatrix);
    if (solution.info() != Success || !solution.eigenvalues().allFinite()){
        throw std::runtime_error("schmidtSplit requires a finite matrix");
    }

    // eigenvalues are in ascending order, and rounding can leave them slightly negative
    VectorXd weights = solution.eigenvalues().cwiseMax(0.0);
    Index size = weights.size();
    double largest = weights(size-1);
    if (largest == 0.0){
        throw std::runtime_error("schmidtSplit requires a nonzero matrix");
    }

    double totalWeight = weights.sum();
    double discarded = 0.0;
    Index keep = size;
    while (keep > 1){
        double weight = weights(size-keep);
        bool tooMany = keep > maxKeep;
        bool isNoise = weight < NOISE_FLOOR*largest;
        bool withinCutoff = discarded + weight <= cutoff*totalWeight;
        if (!tooMany && !isNoise && !withinCutoff){
            break;
        }
        discarded += weight;
        keep--;
    }

    VectorXd kept = weights.tail(keep).reverse().cwiseSqrt();
    kept.normalize();
    MatrixXcd isometry = solution.eigenvectors().rightCols(keep).rowwise().reverse();

    SchmidtSplit split;
    split.schmidtValues = kept;
    split.discardedWeight = discarded/totalWeight;
    if (centerRight){
        split.left = isometry;
        split.right = isometry.adjoint()*matrix;
        split.right /= split.right.norm();
    } else {
        split.left = matrix*isometry;
        split.left /= split.left.norm();
        split.right = isometry.adjoint();
    }
    return split;
}

} //namespace spinchain
