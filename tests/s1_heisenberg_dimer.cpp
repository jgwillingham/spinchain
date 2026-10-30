#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <ExactDiagonalization.hpp>

using namespace Eigen;
using namespace spinchain;



// Two spin-1 sites with isotropic coupling J have H = J/2 [Stot(Stot+1) - 4],
// so the spectrum is -2 (singlet), -1 (triplet, x3), +1 (quintet, x5)
int test_spin1_dimer()
{
        std::vector<double> couplings = {1.0, 1.0, 1.0};
        ExactDiagonalization myspinchain = ExactDiagonalization(SpinChain(2, couplings, 1.0));
        myspinchain.buildHamiltonian();
        myspinchain.solve();

        VectorXcd myspectrum = myspinchain.spectrum;
        assert( myspinchain.dim == 9 );

        double expected[9] = {-2.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        for (int k=0; k<9; k++){
                assert( std::abs( myspectrum[k] - expected[k] ) < 1e-12 );
        }
        return 0;
}



int main(){
	test_spin1_dimer();
}
