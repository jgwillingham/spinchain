#pragma once

#include <eigen3/Eigen/Dense>
#include <Solver.hpp>

namespace spinchain{

using namespace Eigen;


/*! \brief Exact diagonalization of a SpinChain.
 *
 *  Builds the dense many-body Hamiltonian and computes its full spectrum and
 *  eigenstates. The Hilbert space dimension is (2S+1)^L, so chains must be short.
 */
class ExactDiagonalization : public Solver {
    public:
        explicit ExactDiagonalization(const SpinChain& model);

        unsigned int dim;               // Hilbert space dimension (2S+1)^chainlength
        MatrixXcd hamiltonian;          // Hamiltonian
        VectorXcd spectrum;             // Eigenvalues of the Hamiltonian
        MatrixXcd eigenstates;          // Eigenvectors of the Hamiltonian

        // construct the many-body Hamiltonian
        void buildHamiltonian();

        // solve the eigensystem, building the Hamiltonian first if needed
        SolveResult solve() override;

        std::vector<std::complex<double>> localExpectations(const MatrixXcd& op) const override;

    private:
        // embeds an operator acting on sites [i, i+width) into the full Hilbert space
        MatrixXcd embed(const MatrixXcd& op, size_t i, size_t width) const;
};


} // namespace spinchain
