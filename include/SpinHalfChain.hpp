#include <vector>
#include <eigen3/Eigen/Dense>
#include <SpinHalf.hpp>

namespace SpinChain{

using namespace Eigen;


/*! \brief Class for exact diagonalization of 1D spin 1/2 models. 
 *  
 *  SpinHalfChain is for building Hamiltonians for spin 1/2 chains and then
 *  computing their spectrum and eigenstates using exact diagonalization
 */
class SpinHalfChain {
    public:
        SpinHalfChain(size_t chainlength, std::vector<double> couplings);

        unsigned int dim;               // Hilbert space dimension
        std::vector<double> couplings;   // Exchange couplings: Jx, Jy, Jz
        double Jx, Jy, Jz;
        std::vector<SpinHalf> chain;    // Chain of spin-1/2 objects
        MatrixXcd hamiltonian;          // Hamiltonian
        VectorXcd spectrum;             // Eigenvalues of the Hamiltonian
        MatrixXcd eigenstates;          // Eigenvectors of the Hamiltonian

        // construct the many-body Hamiltonian
        void buildHamiltonian();

        // solve the eigensystem
        void diagonalize();

        // returns the expected value of the sum of spins as a 3D vector
        Vector3cd getMagnetization(VectorXcd state);

};


} // namespace SpinChain
