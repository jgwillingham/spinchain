#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <SpinChain.hpp>

using namespace Eigen;
using namespace spinchain;


// decoupled spins in a field h along `axis` have ground energy -L*S*|h|
// and the ground state is fully polarized along the field
int test_decoupled(size_t chainlength, double spin, std::vector<double> field, std::vector<double> expectedDirection)
{
        std::vector<double> couplings = {0.0, 0.0, 0.0};
        SpinChain chain = SpinChain(chainlength, couplings, spin, field);
        chain.buildHamiltonian();
        chain.diagonalize();

        double hmag = std::sqrt(field[0]*field[0] + field[1]*field[1] + field[2]*field[2]);
        assert( std::abs(chain.spectrum[0].real() - (-1.0*chainlength*spin*hmag)) < 1e-12 );

        VectorXcd gs = chain.eigenstates.col(0);
        Vector3cd mag = chain.getMagnetization(gs);
        for (int k=0; k<3; k++){
                assert( std::abs(mag[k].real() - chainlength*spin*expectedDirection[k]) < 1e-12 );
        }
        return 0;
}


// zero field must leave the Hamiltonian unchanged
int test_zero_field()
{
        std::vector<double> couplings = {1.0, 1.0, 1.0};
        SpinChain a = SpinChain(3, couplings, 1.0);
        SpinChain b = SpinChain(3, couplings, 1.0, {0.0, 0.0, 0.0});
        a.buildHamiltonian();
        b.buildHamiltonian();
        assert( (a.hamiltonian - b.hamiltonian).norm() < 1e-15 );
        return 0;
}


int main(){
        test_decoupled(3, 0.5, {0.0, 0.0, 2.0}, {0.0, 0.0, 1.0});
        test_decoupled(2, 1.0, {3.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
        test_decoupled(2, 1.5, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.0});
        test_decoupled(2, 0.5, {1.0, 1.0, 1.0}, {1/std::sqrt(3.0), 1/std::sqrt(3.0), 1/std::sqrt(3.0)});
        test_zero_field();
}
