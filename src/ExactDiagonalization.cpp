#include <ExactDiagonalization.hpp>
#include <math.h>
#include <stdexcept>
#include <eigen3/unsupported/Eigen/KroneckerProduct>


namespace spinchain{


ExactDiagonalization::ExactDiagonalization(const SpinChain& model) : Solver(model)
{
    // dimension of the many body Hilbert space
    dim = pow(model.site.dim, model.chainlength);
}



MatrixXcd ExactDiagonalization::embed(const MatrixXcd& op, size_t i, size_t width) const
{
    size_t dimLeft = pow(model.site.dim, i);
    size_t dimRight = pow(model.site.dim, model.chainlength - i - width);

    MatrixXcd left = kroneckerProduct(MatrixXcd::Identity(dimLeft,dimLeft), op);
    return kroneckerProduct(left, MatrixXcd::Identity(dimRight,dimRight));
}



void ExactDiagonalization::buildHamiltonian()
{
    hamiltonian = MatrixXcd::Zero(dim,dim);

    // the two-site coupling is the same on every bond
    unsigned int bondDim = model.site.dim*model.site.dim;
    MatrixXcd bond = MatrixXcd::Zero(bondDim,bondDim);
    for (const BondTerm& term : model.bondTerms()){
        bond += term.coupling*kroneckerProduct(term.left, term.right);
    }

    for (size_t i=0; i+1<model.chainlength; i++){
        hamiltonian += embed(bond, i, 2);
    }

    // Zeeman coupling to the external field on every site
    MatrixXcd zeeman = model.onsiteTerm();
    for (size_t i=0; i<model.chainlength; i++){
        hamiltonian += embed(zeeman, i, 1);
    }
}



SolveResult ExactDiagonalization::solve()
{
    if (hamiltonian.size() == 0){
        buildHamiltonian();
    }

    SelfAdjointEigenSolver<MatrixXcd> solution(hamiltonian);
    spectrum = solution.eigenvalues();
    eigenstates = solution.eigenvectors();

    return {spectrum[0].real(), true, 1};
}



std::vector<std::complex<double>> ExactDiagonalization::localExpectations(const MatrixXcd& op) const
{
    if (eigenstates.size() == 0){
        throw std::logic_error("localExpectations requires solve() to be called first");
    }

    VectorXcd groundstate = eigenstates.col(0);
    std::vector<std::complex<double>> values;
    for (size_t i=0; i<model.chainlength; i++){
        values.push_back(groundstate.dot(embed(op, i, 1)*groundstate));
    }
    return values;
}

} //namespace spinchain
