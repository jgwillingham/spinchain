#include <SpinChain.hpp>
#include <math.h>
#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace spinchain{


/*! \brief Constructor for a chain of spin-S objects
 *
 *  The spin magnitude is shared by every site, so a single Spin object holds
 *  the local operators for the whole chain.
 */
SpinChain::SpinChain(size_t chainlength, std::vector<double> couplings, double spin,
                     std::vector<double> field)
    : chainlength(chainlength), site(spin)
{
    // save couplings
    this->couplings = couplings;
    Jx = couplings[0];
    Jy = couplings[1];
    Jz = couplings[2];

    // save external field
    this->field = field;
    hx = field[0];
    hy = field[1];
    hz = field[2];

    // dimension of the many body Hilbert space
    dim = pow(site.dim, chainlength);
}



MatrixXcd SpinChain::embed(const MatrixXcd& op, size_t i, size_t width)
{
    size_t dimLeft = pow(site.dim, i);
    size_t dimRight = pow(site.dim, chainlength - i - width);

    MatrixXcd left = kroneckerProduct(MatrixXcd::Identity(dimLeft,dimLeft), op);
    return kroneckerProduct(left, MatrixXcd::Identity(dimRight,dimRight));
}



void SpinChain::buildHamiltonian()
{
    hamiltonian = MatrixXcd::Zero(dim,dim);

    // the two-site coupling is the same on every bond
    MatrixXcd bond = Jx*kroneckerProduct(site.Sx, site.Sx)
                   + Jy*kroneckerProduct(site.Sy, site.Sy)
                   + Jz*kroneckerProduct(site.Sz, site.Sz);

    for (size_t i=0; i+1<chainlength; i++){
        hamiltonian += embed(bond, i, 2);
    }

    // Zeeman coupling to the external field on every site
    MatrixXcd zeeman = -(hx*site.Sx + hy*site.Sy + hz*site.Sz);
    for (size_t i=0; i<chainlength; i++){
        hamiltonian += embed(zeeman, i, 1);
    }
}



void SpinChain::diagonalize()
{
    SelfAdjointEigenSolver<MatrixXcd> solution(hamiltonian);
    spectrum = solution.eigenvalues();
    eigenstates = solution.eigenvectors();
}



Vector3cd SpinChain::getMagnetization(VectorXcd state)
{
    MatrixXcd Mx = MatrixXcd::Zero(dim,dim);
    MatrixXcd My = MatrixXcd::Zero(dim,dim);
    MatrixXcd Mz = MatrixXcd::Zero(dim,dim);

    for (size_t i=0; i<chainlength; i++){
        Mx += embed(site.Sx, i, 1);
        My += embed(site.Sy, i, 1);
        Mz += embed(site.Sz, i, 1);
    }
    std::complex<double> expectedMx, expectedMy, expectedMz;
    expectedMx = state.conjugate().transpose() * Mx * state;
    expectedMy = state.conjugate().transpose() * My * state;
    expectedMz = state.conjugate().transpose() * Mz * state;

    Vector3cd magnetization;
    magnetization << expectedMx, expectedMy, expectedMz;
    return magnetization;
}

} //namespace spinchain
