#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinHalfChain.hpp>

using namespace Eigen;
using namespace SpinChain;


int test_y_hamiltonian()
{
        std::vector<double> couplings = {0.0, 1.0, 0.0};
        SpinHalfChain shc = SpinHalfChain(2, couplings);
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
