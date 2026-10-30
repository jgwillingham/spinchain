#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>
#include <MPO.hpp>
#include <MPS.hpp>

// Expands MPS and MPO objects to dense form, without using any of the library's
// contractions, so tests of those contractions have something independent to compare to.
// Basis ordering matches ExactDiagonalization: site 0 is the most significant digit.

using namespace Eigen;
using namespace spinchain;


// the digits of a basis state index, one for each site
static std::vector<Index> basisDigits(Index state, Index d, size_t length)
{
        std::vector<Index> digits(length);
        for (size_t i=length; i>0; i--){
                digits[i-1] = state % d;
                state /= d;
        }
        return digits;
}


static VectorXcd denseState(const MPS& mps)
{
        Index d = mps.tensor(0).dimension(1);
        Index dim = 1;
        for (size_t i=0; i<mps.length(); i++){
                dim *= d;
        }

        VectorXcd state(dim);
        for (Index k=0; k<dim; k++){
                std::vector<Index> s = basisDigits(k, d, mps.length());
                MatrixXcd product = MatrixXcd::Identity(1,1);
                for (size_t i=0; i<mps.length(); i++){
                        const Tensor3& A = mps.tensor(i);
                        MatrixXcd slice(A.dimension(0), A.dimension(2));
                        for (Index l=0; l<slice.rows(); l++){
                                for (Index r=0; r<slice.cols(); r++){
                                        slice(l,r) = A(l, s[i], r);
                                }
                        }
                        product = (product*slice).eval();
                }
                state(k) = product(0,0);
        }
        return state;
}


static MatrixXcd denseOperator(const MPO& mpo)
{
        size_t length = mpo.W.size();
        Index d = mpo.W[0].dimension(2);
        Index dim = 1;
        for (size_t i=0; i<length; i++){
                dim *= d;
        }

        MatrixXcd op(dim,dim);
        for (Index row=0; row<dim; row++){
                std::vector<Index> s = basisDigits(row, d, length);
                for (Index col=0; col<dim; col++){
                        std::vector<Index> t = basisDigits(col, d, length);
                        MatrixXcd product = MatrixXcd::Identity(1,1);
                        for (size_t i=0; i<length; i++){
                                const Tensor4& W = mpo.W[i];
                                MatrixXcd slice(W.dimension(0), W.dimension(1));
                                for (Index a=0; a<slice.rows(); a++){
                                        for (Index b=0; b<slice.cols(); b++){
                                                slice(a,b) = W(a, b, s[i], t[i]);
                                        }
                                }
                                product = (product*slice).eval();
                        }
                        op(row,col) = product(0,0);
                }
        }
        return op;
}
