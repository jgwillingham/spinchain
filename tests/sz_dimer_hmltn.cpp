#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinHalfChain.hpp>

using namespace Eigen;
using namespace SpinChain;


int test_z_hamiltonian()
{
        std::vector<double> couplings = {0.0, 0.0, 1.0};
        SpinHalfChain shc = SpinHalfChain(2, couplings);
        shc.buildHamiltonian();

        MatrixXcd zHamiltonian;
        zHamiltonian.resize(4, 4);
        zHamiltonian << 0.25, 0.0, 0.0, 0.0,
                        0.0, -0.25, 0.0, 0.0,
                        0.0, 0.0, -0.25, 0.0,
                        0.0, 0.0, 0.0, 0.25;

        assert( shc.hamiltonian.isApprox( zHamiltonian ) );
        return 0;
}


int main(){
	test_z_hamiltonian();
}
