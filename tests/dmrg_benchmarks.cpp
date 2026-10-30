#include <iostream>
#include <assert.h>
#include <cmath>
#include <DMRG.hpp>

using namespace Eigen;
using namespace spinchain;

static const double PI = 3.14159265358979323846;


static DMRGOptions withChi(Index chi)
{
        DMRGOptions options;
        options.chiMax = chi;
        options.energyTol = 1e-6;
        return options;
}


// The open spin-1/2 XX chain in a field maps to free fermions, with single particle
// energies e_k = J cos(pi k/(L+1)) - hz. The ground state fills every negative level.
static double exactXXEnergy(size_t chainlength, double J, double hz)
{
        double total = chainlength*hz/2.0;
        for (size_t k=1; k<=chainlength; k++){
                total += std::min(0.0, J*std::cos(PI*k/(chainlength + 1.0)) - hz);
        }
        return total;
}


// This chain is far too long for exact diagonalization, and critical, so it is a hard
// test. A field also moves the ground state out of the Sz=0 sector, which DMRG could
// never reach from a starting state with a definite Sz.
int test_xx_chain(size_t chainlength, double hz, Index chi, double tolerance)
{
        SpinChain model(chainlength, {1.0, 1.0, 0.0}, 0.5, {0.0, 0.0, hz});
        DMRG dmrg(model, withChi(chi));
        SolveResult result = dmrg.solve();
        double exact = exactXXEnergy(chainlength, 1.0, hz);

        // variational, up to the accuracy of the contraction
        assert( result.energy > exact - 1e-9 );
        assert( result.energy - exact < tolerance );
        return 0;
}


// Fully polarized by a strong field: all the weight must leave the Sz=0 sector
int test_polarized_by_field()
{
        size_t chainlength = 40;
        double hz = 2.5;
        SpinChain model(chainlength, {1.0, 1.0, 1.0}, 0.5, {0.0, 0.0, hz});
        DMRG dmrg(model, withChi(32));
        SolveResult result = dmrg.solve();

        // the polarized state is an exact eigenstate of the Heisenberg chain
        double exact = (chainlength - 1)*0.25 - chainlength*hz/2.0;
        assert( std::abs(result.energy - exact) < 1e-8 );
        assert( std::abs(dmrg.magnetization()[2].real() - chainlength/2.0) < 1e-6 );
        return 0;
}


// The Ising ferromagnet has a diagonal Hamiltonian, so every product state in the z-basis
// is an eigenstate. DMRG must still find the true ground state, whatever the seed.
int test_ising_ferromagnet()
{
        size_t chainlength = 30;
        for (unsigned int seed : {1, 2, 3, 4, 5}){
                DMRGOptions options = withChi(16);
                options.seed = seed;
                DMRG dmrg(SpinChain(chainlength, {0.0, 0.0, -1.0}, 0.5), options);
                assert( std::abs(dmrg.solve().energy - (-0.25*(chainlength - 1))) < 1e-8 );
        }
        return 0;
}


// The bulk energy per site comes from comparing two lengths, which cancels the edge energy.
// The two values are the Bethe ansatz result 1/4 - ln 2, and the numerical result for S=1.
int test_bulk_energy(double spin, size_t shorterLength, size_t longerLength, Index chi, double expected, double tolerance)
{
        std::vector<double> heisenberg = {1.0, 1.0, 1.0};
        double shorter = DMRG(SpinChain(shorterLength, heisenberg, spin), withChi(chi)).solve().energy;
        double longer = DMRG(SpinChain(longerLength, heisenberg, spin), withChi(chi)).solve().energy;
        assert( std::abs((longer - shorter)/(longerLength - shorterLength) - expected) < tolerance );
        return 0;
}


// The energy must fall as the bond dimension grows, with a matching fall in truncation error
int test_chi_dependence()
{
        SpinChain model(20, {1.0, 1.0, 1.0}, 0.5);
        double previousEnergy = 0.0;
        double previousError = 1.0;
        for (Index chi : {8, 16, 32}){
                DMRG dmrg(model, withChi(chi));
                double result = dmrg.solve().energy;
                assert( result < previousEnergy + 1e-9 );
                assert( dmrg.maxTruncationError() < previousError );
                previousEnergy = result;
                previousError = dmrg.maxTruncationError();
        }
        return 0;
}


// The energy of a gapped chain must not depend on the random starting state
int test_seed_independence()
{
        SpinChain model(16, {1.0, 0.9, 1.3}, 0.5, {0.1, 0.0, 0.2});
        std::vector<double> energies;
        for (unsigned int seed : {1, 2, 3}){
                DMRGOptions options = withChi(32);
                options.seed = seed;
                energies.push_back(DMRG(model, options).solve().energy);
        }
        assert( std::abs(energies[0] - energies[1]) < 1e-6 );
        assert( std::abs(energies[0] - energies[2]) < 1e-6 );
        return 0;
}


int main(){
        test_xx_chain(40, 0.0, 32, 1e-5);
        test_xx_chain(40, 0.4, 32, 1e-5);
        test_xx_chain(40, 1.5, 16, 1e-8);
        test_polarized_by_field();
        test_ising_ferromagnet();
        test_bulk_energy(0.5, 30, 40, 32, 0.25 - std::log(2.0), 1e-3);
        // shorter S=1 chains fail this: their edge spins still interact through the bulk
        test_bulk_energy(1.0, 30, 40, 32, -1.401484, 1e-3);
        test_chi_dependence();
        test_seed_independence();
}
