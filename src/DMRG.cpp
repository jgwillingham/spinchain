#include <DMRG.hpp>
#include <Lanczos.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>


namespace spinchain{

// bond dimension of the random starting state
static const Index INITIAL_BOND_DIM = 8;

// loose two-site solves are enough while the environments are still poor;
// the Lanczos tolerance starts here and tightens by a factor of 10 each sweep
static const double INITIAL_LANCZOS_TOL = 1e-3;


// the bond dimension doubles every sweep until it reaches chiMax
static Index chiForSweep(size_t sweep, Index chiMax)
{
    Index chi = std::min(2*INITIAL_BOND_DIM, chiMax);
    for (size_t k=0; k<sweep && chi<chiMax; k++){
        chi = std::min(2*chi, chiMax);
    }
    return chi;
}


static double lanczosTolForSweep(size_t sweep, double finalTol)
{
    return std::max(finalTol, INITIAL_LANCZOS_TOL*std::pow(0.1, static_cast<double>(sweep)));
}


static DMRGOptions validated(const SpinChain& model, const DMRGOptions& options)
{
    if (model.chainlength < 2){
        throw std::invalid_argument("DMRG requires at least 2 sites");
    }
    if (options.chiMax < 1 || options.maxSweeps < 1 || options.lanczosIterations < 1){
        throw std::invalid_argument("chiMax, maxSweeps and lanczosIterations must be positive");
    }
    return options;
}



DMRG::DMRG(const SpinChain& model, DMRGOptions options)
    : Solver(model),
      options(validated(model, options)),
      hamiltonian(model),
      state(model.chainlength, model.site.dim, std::min(INITIAL_BOND_DIM, options.chiMax), options.seed),
      left(model.chainlength),
      right(model.chainlength),
      truncationError(0.0),
      solved(false)
{
}



void DMRG::optimizeBond(size_t i, Index chi, double lanczosTol, bool centerRight)
{
    Tensor4 theta = state.twoSiteTensor(i);
    Eigen::array<Index,4> dims = theta.dimensions();

    auto applyH = [&](const VectorXcd& x) -> VectorXcd {
        Tensor4 result = applyTwoSiteHamiltonian(left[i], hamiltonian.W[i], hamiltonian.W[i+1],
                                                 right[i+1], TensorMap<const Tensor4>(x.data(), dims));
        return Map<VectorXcd>(result.data(), result.size());
    };

    // the current state is an excellent starting vector, since it is close to converged
    VectorXcd start = Map<VectorXcd>(theta.data(), theta.size());
    LanczosResult lowest = lanczos(applyH, start, options.lanczosIterations, lanczosTol);

    Tensor4 optimized = TensorMap<Tensor4>(lowest.eigenvector.data(), dims);
    double discarded = state.setTwoSiteTensor(i, optimized, chi, options.truncCutoff, centerRight);
    truncationError = std::max(truncationError, discarded);
}



SolveResult DMRG::solve()
{
    size_t length = model.chainlength;

    state = MPS(length, model.site.dim, std::min(INITIAL_BOND_DIM, options.chiMax), options.seed);

    left[0] = emptyEnvironment();
    right[length-1] = emptyEnvironment();
    for (size_t i=length-1; i>0; i--){
        right[i-1] = growRight(right[i], state.tensor(i), hamiltonian.W[i]);
    }

    double energyNow = std::numeric_limits<double>::infinity();
    bool converged = false;
    size_t sweep = 0;

    while (sweep < options.maxSweeps && !converged){
        Index chi = chiForSweep(sweep, options.chiMax);
        double lanczosTol = lanczosTolForSweep(sweep, options.lanczosTol);
        truncationError = 0.0;

        for (size_t i=0; i+1<length; i++){
            optimizeBond(i, chi, lanczosTol, true);
            left[i+1] = growLeft(left[i], state.tensor(i), hamiltonian.W[i]);
        }
        for (size_t i=length-1; i>0; i--){
            optimizeBond(i-1, chi, lanczosTol, false);
            right[i-1] = growRight(right[i], state.tensor(i), hamiltonian.W[i]);
        }

        double energyBefore = energyNow;
        energyNow = energy(state, hamiltonian);

        // an energy change only means something between sweeps with the final settings.
        // The truncation error plays no part here: for a critical chain it stays above
        // the cutoff at any chiMax, and it is reported to the caller instead.
        bool finalSettings = sweep > 0 && chi == options.chiMax
                          && chiForSweep(sweep-1, options.chiMax) == chi
                          && lanczosTol == options.lanczosTol;
        converged = finalSettings && std::abs(energyNow - energyBefore) < options.energyTol;

        sweep++;
        if (options.verbose){
            std::cout << "sweep " << sweep << ": chi = " << state.maxBondDim()
                      << ", E = " << energyNow
                      << ", truncation error = " << truncationError << std::endl;
        }
    }

    solved = true;
    return {energyNow, converged, sweep};
}



std::vector<std::complex<double>> DMRG::localExpectations(const MatrixXcd& op) const
{
    if (!solved){
        throw std::logic_error("localExpectations requires solve() to be called first");
    }
    return state.localExpectations(op);
}



double DMRG::entanglementEntropy(size_t bond) const
{
    if (!solved){
        throw std::logic_error("entanglementEntropy requires solve() to be called first");
    }
    return state.entanglementEntropy(bond);
}



double DMRG::maxTruncationError() const
{
    return truncationError;
}



const MPS& DMRG::groundState() const
{
    return state;
}

} //namespace spinchain
