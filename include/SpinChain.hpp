#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>
#include <Spin.hpp>

namespace spinchain{

using namespace Eigen;


/*! \brief Class for exact diagonalization of 1D spin-S models.
 *
 *  SpinChain is for building Hamiltonians for open spin-S chains with
 *  nearest-neighbor XYZ couplings and then computing their spectrum and
 *  eigenstates using exact diagonalization.
 */
class SpinChain {
    public:
        /*! field is the external magnetic field (hx, hy, hz); it enters the
         *  Hamiltonian as the Zeeman term -sum_i (hx Sx_i + hy Sy_i + hz Sz_i).
         */
        SpinChain(size_t chainlength, std::vector<double> couplings, double spin,
                  std::vector<double> field = {0.0, 0.0, 0.0});

        size_t chainlength;             // Number of sites
        Spin site;                      // Single-site spin operators
        unsigned int dim;               // Hilbert space dimension (2S+1)^chainlength
        std::vector<double> couplings;  // Exchange couplings: Jx, Jy, Jz
        double Jx, Jy, Jz;
        std::vector<double> field;      // External magnetic field: hx, hy, hz
        double hx, hy, hz;
        MatrixXcd hamiltonian;          // Hamiltonian
        VectorXcd spectrum;             // Eigenvalues of the Hamiltonian
        MatrixXcd eigenstates;          // Eigenvectors of the Hamiltonian

        // construct the many-body Hamiltonian
        void buildHamiltonian();

        // solve the eigensystem
        void diagonalize();

        // returns the expected value of the sum of spins as a 3D vector
        Vector3cd getMagnetization(VectorXcd state);

    private:
        // embeds an operator acting on sites [i, i+width) into the full Hilbert space
        MatrixXcd embed(const MatrixXcd& op, size_t i, size_t width);
};


} // namespace spinchain
