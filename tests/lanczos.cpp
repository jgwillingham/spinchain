#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <Lanczos.hpp>

using namespace Eigen;
using namespace spinchain;


static MatrixXcd randomHermitian(int n)
{
        MatrixXcd M = MatrixXcd::Random(n,n);
        return 0.5*(M + M.adjoint());
}


// checks a Lanczos result against the exact lowest eigenvalue of H
static void check_against_exact(const MatrixXcd& H, const LanczosResult& result)
{
        SelfAdjointEigenSolver<MatrixXcd> exact(H);
        assert( result.converged );
        assert( std::abs(result.eigenvalue - exact.eigenvalues()(0)) < 1e-9 );
        assert( std::abs(result.eigenvector.norm() - 1.0) < 1e-12 );
        assert( (H*result.eigenvector - result.eigenvalue*result.eigenvector).norm() < 1e-8 );
}


// a small matrix is solved exactly once the Krylov space fills the whole space
int test_random_matrix()
{
        for (int n : {1, 2, 3, 8, 40}){
                MatrixXcd H = randomHermitian(n);
                auto applyH = [&H](const VectorXcd& x) -> VectorXcd { return H*x; };
                check_against_exact(H, lanczos(applyH, VectorXcd::Random(n), n, 1e-10));
        }
        return 0;
}


// degenerate lowest level: any vector in the eigenspace is a correct answer
int test_degenerate()
{
        VectorXd levels(6);
        levels << -2.0, -2.0, -2.0, 0.5, 1.0, 3.0;
        MatrixXcd Q = HouseholderQR<MatrixXcd>(MatrixXcd::Random(6,6)).householderQ();
        MatrixXcd H = Q*levels.cast<std::complex<double>>().asDiagonal()*Q.adjoint();
        auto applyH = [&H](const VectorXcd& x) -> VectorXcd { return H*x; };
        check_against_exact(H, lanczos(applyH, VectorXcd::Random(6), 6, 1e-10));
        return 0;
}


// starting on an exact eigenvector must stop immediately, without dividing by a zero beta
int test_breakdown()
{
        MatrixXcd H = MatrixXcd::Zero(4,4);
        H.diagonal() << 1.0, 2.0, 3.0, 4.0;
        auto applyH = [&H](const VectorXcd& x) -> VectorXcd { return H*x; };

        VectorXcd start = VectorXcd::Unit(4, 2);
        LanczosResult result = lanczos(applyH, start, 4, 0.0);
        assert( result.converged );
        assert( result.iterations == 1 );
        assert( std::abs(result.eigenvalue - 3.0) < 1e-14 );
        return 0;
}


// running out of iterations must be reported rather than hidden
int test_not_converged()
{
        MatrixXcd H = randomHermitian(30);
        auto applyH = [&H](const VectorXcd& x) -> VectorXcd { return H*x; };
        LanczosResult result = lanczos(applyH, VectorXcd::Random(30), 3, 1e-12);
        assert( !result.converged );
        assert( result.iterations == 3 );
        return 0;
}


int main(){
        test_random_matrix();
        test_degenerate();
        test_breakdown();
        test_not_converged();
}
