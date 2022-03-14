#include <SpinHalfChain.hpp>
#include <math.h>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <iostream>


namespace SpinChain{

using namespace std::literals;


/*! \brief Constructor for a chain of spin-1/2 objects
 *  
 *  When a SpinHalfChain object is declared, a vector of SpinHalf objects is assembled
 *  for the given chain length.
 */
SpinHalfChain::SpinHalfChain(size_t chainlength, std::vector<float> couplings)
{
    // save couplings
    couplings = couplings;

    // assemble a vector of spin-1/2 objects
    chain.assign(chainlength, SpinHalf());

    dim = pow(2, chainlength); // dimension of the many body Hilbert space
}




void SpinHalfChain::buildHamiltonian()
{
    hamiltonian.resize(dim,dim);

    size_t chainlength = chain.size();

    for (int i=0; i<chainlength-1; i++){
        MatrixXcf SxSx, SySy, SzSz;
        // the i and i+1 couplings
        SxSx = kroneckerProduct(chain[i].Sx, chain[i].Sx);
        SySy = kroneckerProduct(chain[i].Sy, chain[i].Sy);
        SzSz = kroneckerProduct(chain[i].Sz, chain[i].Sz);

        for (int j=i-1; j>=0; j--){
            SxSx = kroneckerProduct(Matrix2cf::Identity(), SxSx).eval();
            SySy = kroneckerProduct(Matrix2cf::Identity(), SySy).eval();
            SzSz = kroneckerProduct(Matrix2cf::Identity(), SzSz).eval();
        }
        for (int j=i+2; j<chainlength; j++){
            SxSx = kroneckerProduct(SxSx, Matrix2cf::Identity()).eval();
            SySy = kroneckerProduct(SySy, Matrix2cf::Identity()).eval();
            SzSz = kroneckerProduct(SzSz, Matrix2cf::Identity()).eval();    
        }

        hamiltonian += (SxSx + SySy + SzSz);
    }
}



void SpinHalfChain::diagonalize()
{
    SelfAdjointEigenSolver<MatrixXcf> solution(hamiltonian);
    spectrum = solution.eigenvalues();
    eigenstates = solution.eigenvectors();
}

} //namespace SpinChain