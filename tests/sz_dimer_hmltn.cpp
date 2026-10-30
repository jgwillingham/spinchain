#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinChain.hpp>

using namespace Eigen;
using namespace spinchain;


int test_z_hamiltonian()
{
        std::vector<double> couplings = {0.0, 0.0, 1.0};
        SpinChain shc = SpinChain(2, couplings, 0.5);
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
