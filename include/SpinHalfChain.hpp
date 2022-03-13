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
        SpinHalfChain(unsigned int length, std::vector<float> couplings);

        // Exchange couplings (Jx,Jy,Jz) for the model
        std::vector<float> couplings;

        // Chain of spin-1/2 objects
        std::vector<SpinHalf> chain;

        // Hamiltonian
        MatrixXf hamiltonian;

        // eigenvalues of the model Hamiltonian
        VectorXf spectrum;

        // eigenvectors of the model Hamiltonian
        MatrixXf eigenstates;

        // construct the many-body Hamiltonian
        void buildHamiltonian();

        // solve the eigensystem
        void diagonalize();

        // returns the expected value of the sum of spins as a 3D vector
        Vector3f getMagnetization();

};


} // namespace SpinChain