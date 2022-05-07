#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinHalfChain.hpp>

using namespace Eigen;
using namespace SpinChain;



int test_groundstate(size_t chainlength, std::vector<double> couplings)
{

        SpinHalfChain myspinchain = SpinHalfChain(chainlength, couplings);
        myspinchain.buildHamiltonian();
        myspinchain.diagonalize();

        VectorXcd myspectrum;
        myspectrum = myspinchain.spectrum;

        //std::cout << myspectrum << std::endl;
        assert( myspectrum[0].imag() == 0.0 );
        assert( std::abs( myspectrum[0].real() - -0.25 ) < 1e-15 );
        assert( myspectrum[0].real() == (-1.0*myspectrum[3].real()) );
        return 0;
}




int main(){
	size_t chainlength = 2;
        std::vector<double> couplings = {0.0, 0.0, 1.0};
	test_groundstate(chainlength, couplings);
}
