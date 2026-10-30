#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinChain.hpp>

using namespace Eigen;
using namespace spinchain;



int test_x_hamiltonian()
{
	std::vector<double> couplings = {1.0, 0.0, 0.0};
	SpinChain shc = SpinChain(2, couplings, 0.5);
	shc.buildHamiltonian();

	MatrixXcd xHamiltonian;
	xHamiltonian.resize(4, 4);
	xHamiltonian << 0.0, 0.0, 0.0, 0.25,
		        0.0, 0.0, 0.25, 0.0,
		        0.0, 0.25, 0.0, 0.0,
			0.25, 0.0, 0.0, 0.0;	

	assert( shc.hamiltonian.isApprox( xHamiltonian ) );
	return 0;
}


int main(){
	test_x_hamiltonian();
}
