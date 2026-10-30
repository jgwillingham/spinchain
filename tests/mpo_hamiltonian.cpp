#include <iostream>
#include <assert.h>
#include <random>
#include <ExactDiagonalization.hpp>
#include <MPO.hpp>
#include "dense_helpers.hpp"


// the MPO must reproduce the exact Hamiltonian entry by entry
int test_matches_exact(size_t chainlength, double spin, std::vector<double> couplings, std::vector<double> field)
{
        SpinChain model(chainlength, couplings, spin, field);
        ExactDiagonalization exact(model);
        exact.buildHamiltonian();

        MatrixXcd fromMPO = denseOperator(MPO(model));
        assert( (fromMPO - exact.hamiltonian).norm() < 1e-12 );
        return 0;
}


int main(){
        std::mt19937 generator(7);
        std::uniform_real_distribution<double> random(-1.0, 1.0);

        for (size_t chainlength=2; chainlength<=6; chainlength++){
                for (double spin : {0.5, 1.0}){
                        std::vector<double> couplings = {random(generator), random(generator), random(generator)};
                        std::vector<double> field = {random(generator), random(generator), random(generator)};
                        test_matches_exact(chainlength, spin, couplings, field);
                }
        }
        test_matches_exact(4, 1.5, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0});
}
