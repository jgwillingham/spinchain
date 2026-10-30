#pragma once

#include <vector>
#include <Environments.hpp>
#include <MPO.hpp>
#include <MPS.hpp>
#include <Solver.hpp>

namespace spinchain{


/*! \brief Settings for a DMRG calculation
 */
struct DMRGOptions {
    Index chiMax = 100;             // largest bond dimension; this sets the accuracy
    size_t maxSweeps = 30;          // a sweep is one pass right and one pass back
    double energyTol = 1e-8;        // stop when the energy changes by less than this in a sweep
    double truncCutoff = 1e-12;     // largest discarded weight at each truncation
    size_t lanczosIterations = 8;   // largest Krylov space for each two-site update; too small
                                    // and observables converge much more slowly than the energy
    double lanczosTol = 1e-8;       // final residual tolerance for each two-site update
    unsigned int seed = 1;          // seed for the random starting state
    bool verbose = false;           // print progress after each sweep
};


/*! \brief Two-site finite DMRG for the ground state of an open SpinChain
 *
 *  DMRG represents the ground state as a matrix product state and improves it by
 *  sweeping along the chain. At each bond, the two sites are solved exactly in the
 *  basis given by the rest of the chain, and then split back apart keeping only the
 *  largest Schmidt values. Optimizing two sites at once lets the bond dimension adapt
 *  to the state.
 *
 *  The bond dimension doubles every sweep up to chiMax. Accuracy is controlled by
 *  chiMax: if it is too small, the calculation still converges, but to an energy
 *  above the true ground energy. Check this by comparing results for several chiMax,
 *  and with maxTruncationError().
 *
 *  If the ground state is degenerate, the result can be any state in that manifold
 *  and depends on the seed. Only the energy is meaningful then.
 *
 *  The tensors are dense, so no symmetry of the model is required or exploited.
 */
class DMRG : public Solver {
    public:
        explicit DMRG(const SpinChain& model, DMRGOptions options = DMRGOptions());

        // returns <psi|H|psi> of the final state, which is variational unlike the
        // eigenvalues from the two-site updates, since those are before truncation
        SolveResult solve() override;

        std::vector<std::complex<double>> localExpectations(const MatrixXcd& op) const override;

        // returns the von Neumann entropy across the bond between sites i and i+1
        double entanglementEntropy(size_t bond) const;

        // returns the largest discarded weight of any truncation in the last sweep
        double maxTruncationError() const;

        const MPS& groundState() const;

    private:
        // solves sites i and i+1 together and leaves the center on the side given by centerRight
        void optimizeBond(size_t i, Index chi, double lanczosTol, bool centerRight);

        DMRGOptions options;
        MPO hamiltonian;
        MPS state;
        std::vector<Tensor3> left;      // left[i] covers all sites before i
        std::vector<Tensor3> right;     // right[i] covers all sites after i
        double truncationError;
        bool solved;
};

} // namespace spinchain
