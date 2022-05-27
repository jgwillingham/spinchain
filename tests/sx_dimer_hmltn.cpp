#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinHalfChain.hpp>

using namespace Eigen;
using namespace SpinChain;



int test_x_hamiltonian()
{
	std::vector<double> couplings = {1.0, 0.0, 0.0};
	SpinHalfChain shc = SpinHalfChain(2, couplings);
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
