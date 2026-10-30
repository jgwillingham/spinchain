#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <eigen3/Eigen/Dense>
#include <SchmidtSplit.hpp>

using namespace Eigen;
using namespace spinchain;


// builds a matrix with the given Schmidt values
static MatrixXcd withSchmidtValues(const VectorXd& s, int rows, int cols)
{
        MatrixXcd U = HouseholderQR<MatrixXcd>(MatrixXcd::Random(rows,rows)).householderQ();
        MatrixXcd V = HouseholderQR<MatrixXcd>(MatrixXcd::Random(cols,cols)).householderQ();
        return U.leftCols(s.size())*s.cast<std::complex<double>>().asDiagonal()*V.leftCols(s.size()).adjoint();
}


// keeping everything reproduces the matrix, whichever side holds the Schmidt values
int test_exact()
{
        MatrixXcd M = MatrixXcd::Random(6,4);
        M /= M.norm();
        for (bool centerRight : {true, false}){
                SchmidtSplit split = schmidtSplit(M, 100, 0.0, centerRight);
                assert( split.schmidtValues.size() == 4 );
                assert( (split.left*split.right - M).norm() < 1e-12 );
                assert( split.discardedWeight < 1e-14 );
        }
        return 0;
}


// the side without the Schmidt values is an isometry and the discarded weight is the sum of the dropped s^2
int test_max_keep()
{
        VectorXd s(4);
        s << 0.8, 0.5, 0.3, 0.1;
        s.normalize();
        MatrixXcd M = withSchmidtValues(s, 7, 5);

        SchmidtSplit toRight = schmidtSplit(M, 2, 0.0, true);
        assert( toRight.schmidtValues.size() == 2 );
        assert( (toRight.schmidtValues - s.head(2)/s.head(2).norm()).norm() < 1e-13 );
        assert( std::abs(toRight.right.norm() - 1.0) < 1e-13 );
        assert( std::abs(toRight.discardedWeight - s.tail(2).squaredNorm()) < 1e-13 );
        assert( (toRight.left.adjoint()*toRight.left - MatrixXcd::Identity(2,2)).norm() < 1e-12 );

        SchmidtSplit toLeft = schmidtSplit(M, 2, 0.0, false);
        assert( std::abs(toLeft.left.norm() - 1.0) < 1e-13 );
        assert( (toLeft.right*toLeft.right.adjoint() - MatrixXcd::Identity(2,2)).norm() < 1e-12 );
        return 0;
}


// the cutoff bounds the total discarded weight, not each value separately
int test_cutoff()
{
        VectorXd s(4);
        s << 1.0, 1e-2, 1e-3, 1e-3;
        s.normalize();
        MatrixXcd M = withSchmidtValues(s, 6, 6);

        // both 1e-3 values fit together (2e-6 of weight), but not with the 1e-2 value
        assert( schmidtSplit(M, 100, 3e-6, true).schmidtValues.size() == 2 );
        // only one of the 1e-3 values fits
        assert( schmidtSplit(M, 100, 1.5e-6, true).schmidtValues.size() == 3 );
        return 0;
}


// values at the level of rounding noise never survive, even with no cutoff
int test_rank_deficient()
{
        VectorXd s(2);
        s << 0.8, 0.6;
        SchmidtSplit split = schmidtSplit(withSchmidtValues(s, 6, 6), 100, 0.0, true);
        assert( split.schmidtValues.size() == 2 );
        return 0;
}


// The spectra of DMRG wavefunctions fall off quickly until they hit machine precision.
// These are the matrices that made Eigen's BDCSVD read out of bounds.
int test_rapidly_decaying_spectrum()
{
        int n = 64;
        VectorXd s(n);
        for (int i=0; i<n; i++){
                s(i) = std::exp(-1.5*i);
        }
        s.normalize();
        MatrixXcd M = withSchmidtValues(s, n, n);

        for (bool centerRight : {true, false}){
                SchmidtSplit split = schmidtSplit(M, n, 1e-14, centerRight);
                Index keep = split.schmidtValues.size();
                assert( keep < n );
                assert( (split.left*split.right - M).norm() < 1e-6 );
                assert( split.discardedWeight <= 1e-14 );
        }
        return 0;
}


// one value always survives, however aggressive the settings
int test_keeps_one()
{
        SchmidtSplit split = schmidtSplit(MatrixXcd::Random(4,4), 0, 1.0, true);
        assert( split.schmidtValues.size() == 1 );
        assert( std::abs(split.schmidtValues(0) - 1.0) < 1e-14 );
        return 0;
}


int test_rejects_bad_input()
{
        bool threw = false;
        try {
                schmidtSplit(MatrixXcd::Zero(3,3), 2, 0.0, true);
        } catch (const std::runtime_error &){
                threw = true;
        }
        assert( threw );
        return 0;
}


int main(){
        test_exact();
        test_max_keep();
        test_cutoff();
        test_rank_deficient();
        test_rapidly_decaying_spectrum();
        test_keeps_one();
        test_rejects_bad_input();
}
