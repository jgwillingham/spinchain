#include <iostream>
#include <assert.h>
#include <ExactDiagonalization.hpp>
#include <Environments.hpp>
#include <Lanczos.hpp>
#include <MPS.hpp>
#include "dense_helpers.hpp"


// A as a matrix with rows (bondLeft, physical) and columns bondRight
static MatrixXcd leftGrouped(const Tensor3& A)
{
        return Map<const MatrixXcd>(A.data(), A.dimension(0)*A.dimension(1), A.dimension(2));
}

// A as a matrix with rows bondLeft and columns (physical, bondRight)
static MatrixXcd rightGrouped(const Tensor3& A)
{
        return Map<const MatrixXcd>(A.data(), A.dimension(0), A.dimension(1)*A.dimension(2));
}


// every tensor left of the center must be a left isometry and every tensor right of it a right isometry
static void check_canonical(const MPS& mps)
{
        for (size_t i=0; i<mps.center(); i++){
                MatrixXcd M = leftGrouped(mps.tensor(i));
                assert( (M.adjoint()*M - MatrixXcd::Identity(M.cols(),M.cols())).norm() < 1e-12 );
        }
        for (size_t i=mps.center()+1; i<mps.length(); i++){
                MatrixXcd M = rightGrouped(mps.tensor(i));
                assert( (M*M.adjoint() - MatrixXcd::Identity(M.rows(),M.rows())).norm() < 1e-12 );
        }
}


// moving the center changes the tensors but must not change the state
int test_center_moves(size_t length, Index d)
{
        MPS mps(length, d, 8, 11);
        check_canonical(mps);

        VectorXcd original = denseState(mps);
        assert( std::abs(original.norm() - 1.0) < 1e-12 );

        for (size_t site : std::vector<size_t>{length-1, 2, 0, 3, 1}){
                mps.moveCenterTo(site);
                assert( mps.center() == site );
                check_canonical(mps);
                assert( (denseState(mps) - original).norm() < 1e-12 );
        }
        return 0;
}


// bonds near the edges are limited by the Hilbert space on their short side
int test_bond_dimensions()
{
        MPS mps(6, 2, 100, 3);
        std::vector<Index> expected = {2, 4, 8, 4, 2, 1};
        for (size_t i=0; i<6; i++){
                assert( mps.tensor(i).dimension(2) == expected[i] );
        }
        assert( mps.maxBondDim() == 8 );
        return 0;
}


// the same seed gives the same state, and different seeds do not
int test_seed()
{
        VectorXcd a = denseState(MPS(4, 2, 4, 1));
        VectorXcd b = denseState(MPS(4, 2, 4, 1));
        VectorXcd c = denseState(MPS(4, 2, 4, 2));
        assert( (a - b).norm() == 0.0 );
        assert( (a - c).norm() > 1e-3 );
        return 0;
}


// compares localExpectations and energy against the dense state and Hamiltonian
int test_expectations(size_t length, double spin)
{
        SpinChain model(length, {0.7, -0.4, 1.1}, spin, {0.3, -0.2, 0.5});
        ExactDiagonalization exact(model);
        exact.buildHamiltonian();

        MPS mps(length, model.site.dim, 6, 5);
        mps.moveCenterTo(length/2);
        VectorXcd psi = denseState(mps);

        std::complex<double> exactEnergy = psi.dot(exact.hamiltonian*psi);
        assert( std::abs(energy(mps, MPO(model)) - exactEnergy.real()) < 1e-12 );

        // Sx is applied to one site at a time on the dense vector
        Index d = model.site.dim;
        std::vector<std::complex<double>> values = mps.localExpectations(model.site.Sx);
        Index dimRight = psi.size();
        Index dimLeft = 1;
        for (size_t i=0; i<length; i++){
                dimRight /= d;
                VectorXcd acted = VectorXcd::Zero(psi.size());
                for (Index l=0; l<dimLeft; l++){
                        for (Index r=0; r<dimRight; r++){
                                for (Index s=0; s<d; s++){
                                        for (Index t=0; t<d; t++){
                                                acted((l*d + s)*dimRight + r) += model.site.Sx(s,t)*psi((l*d + t)*dimRight + r);
                                        }
                                }
                        }
                }
                assert( std::abs(values[i] - psi.dot(acted)) < 1e-12 );
                dimLeft *= d;
        }
        return 0;
}


// with the center on a bond, <theta|H_eff|theta> is the energy of the whole state,
// and with an exact basis the lowest eigenvalue of H_eff is the exact ground energy
int test_effective_hamiltonian()
{
        size_t length = 4;
        size_t bond = 1;
        SpinChain model(length, {1.0, 0.8, 0.6}, 0.5, {0.2, 0.0, -0.3});
        MPO mpo(model);
        MPS mps(length, model.site.dim, 16, 9);
        mps.moveCenterTo(bond);

        Tensor3 left = emptyEnvironment();
        for (size_t i=0; i<bond; i++){
                left = growLeft(left, mps.tensor(i), mpo.W[i]);
        }
        Tensor3 right = emptyEnvironment();
        for (size_t i=length-1; i>bond+1; i--){
                right = growRight(right, mps.tensor(i), mpo.W[i]);
        }

        Tensor4 theta = mps.twoSiteTensor(bond);
        Eigen::array<Index,4> dims = theta.dimensions();
        auto applyH = [&](const VectorXcd& x) -> VectorXcd {
                Tensor4 result = applyTwoSiteHamiltonian(left, mpo.W[bond], mpo.W[bond+1], right,
                                                         TensorMap<const Tensor4>(x.data(), dims));
                return Map<VectorXcd>(result.data(), result.size());
        };

        VectorXcd flat = Map<VectorXcd>(theta.data(), theta.size());
        std::complex<double> expected = flat.dot(applyH(flat));
        assert( std::abs(expected.real() - energy(mps, mpo)) < 1e-12 );
        assert( std::abs(expected.imag()) < 1e-12 );

        ExactDiagonalization exact(model);
        double groundEnergy = exact.solve().energy;
        LanczosResult lowest = lanczos(applyH, flat, flat.size(), 1e-10);
        assert( std::abs(lowest.eigenvalue - groundEnergy) < 1e-9 );
        return 0;
}


int main(){
        test_center_moves(6, 2);
        test_center_moves(5, 3);
        test_bond_dimensions();
        test_seed();
        test_expectations(5, 0.5);
        test_expectations(4, 1.0);
        test_effective_hamiltonian();
}
