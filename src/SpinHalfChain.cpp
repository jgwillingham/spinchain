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
SpinHalfChain::SpinHalfChain(size_t chainlength, std::vector<double> couplings)
{
    // save couplings
    couplings = couplings;
    Jx = couplings[0];
    Jy = couplings[1];
    Jz = couplings[2];

    // assemble a vector of spin-1/2 objects
    chain.assign(chainlength, SpinHalf());

    // dimension of the many body Hilbert space
    dim = pow(2, chainlength); 
}




void SpinHalfChain::buildHamiltonian()
{
    hamiltonian.resize(dim,dim);

    size_t chainlength = chain.size();

    for (int i=0; i<chainlength-1; i++){
        MatrixXcd SxSx, SySy, SzSz;
        // the i and i+1 couplings
        SxSx = kroneckerProduct(chain[i].Sx, chain[i+1].Sx);
        SySy = kroneckerProduct(chain[i].Sy, chain[i+1].Sy);
        SzSz = kroneckerProduct(chain[i].Sz, chain[i+1].Sz);

        for (int j=i-1; j>=0; j--){
            SxSx = kroneckerProduct(Matrix2cd::Identity(), SxSx).eval();
            SySy = kroneckerProduct(Matrix2cd::Identity(), SySy).eval();
            SzSz = kroneckerProduct(Matrix2cd::Identity(), SzSz).eval();
        }
        for (int j=i+2; j<chainlength; j++){
            SxSx = kroneckerProduct(SxSx, Matrix2cd::Identity()).eval();
            SySy = kroneckerProduct(SySy, Matrix2cd::Identity()).eval();
            SzSz = kroneckerProduct(SzSz, Matrix2cd::Identity()).eval();    
        }

        hamiltonian += (Jx*SxSx + Jy*SySy + Jz*SzSz);
    }
}



void SpinHalfChain::diagonalize()
{
    SelfAdjointEigenSolver<MatrixXcd> solution(hamiltonian);
    spectrum = solution.eigenvalues();
    eigenstates = solution.eigenvectors();
}



Vector3cd SpinHalfChain::getMagnetization(VectorXcd state)
{
    int chainlength = chain.size();

    MatrixXcd Mx, My, Mz;
    Mx.resize(dim,dim);
    My.resize(dim,dim);
    Mz.resize(dim,dim);

    for (int i=0; i<chainlength; i++){
        MatrixXcd Mxi = chain[i].Sx;
        MatrixXcd Myi = chain[i].Sy;
        MatrixXcd Mzi = chain[i].Sz;
        for (int j=i-1; j>=0; j--){
            Mxi = kroneckerProduct(Matrix2cd::Identity(), Mxi).eval();
            Myi = kroneckerProduct(Matrix2cd::Identity(), Myi).eval();
            Mzi = kroneckerProduct(Matrix2cd::Identity(), Mzi).eval();
        }
        for (int j=i+1; j<chainlength; j++){
            Mxi = kroneckerProduct(Mxi, Matrix2cd::Identity()).eval();
            Myi = kroneckerProduct(Myi, Matrix2cd::Identity()).eval();
            Mzi = kroneckerProduct(Mzi, Matrix2cd::Identity()).eval();
        }
        Mx += Mxi;
        My += Myi;
        Mz += Mzi;
    }
    std::complex<double> expectedMx, expectedMy, expectedMz;
    expectedMx = state.conjugate().transpose() * Mx * state;
    expectedMy = state.conjugate().transpose() * My * state;
    expectedMz = state.conjugate().transpose() * Mz * state;
    
    Vector3cd magnetization;
    magnetization << expectedMx, expectedMy, expectedMz;
    return magnetization;
}

} //namespace SpinChain
