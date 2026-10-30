#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinChain.hpp>

using namespace Eigen;
using namespace spinchain;


int test_y_hamiltonian()
{
        std::vector<double> couplings = {0.0, 1.0, 0.0};
        SpinChain shc = SpinChain(2, couplings, 0.5);
        shc.buildHamiltonian();

        MatrixXcd yHamiltonian;
        yHamiltonian.resize(4, 4);
        yHamiltonian << 0.0, 0.0, 0.0, -0.25,
                        0.0, 0.0, 0.25, 0.0,
                        0.0, 0.25, 0.0, 0.0,
                        -0.25, 0.0, 0.0, 0.0;

        assert( shc.hamiltonian.isApprox( yHamiltonian ) );
        return 0;
}



int main(){
	test_y_hamiltonian();
}
