#include <iostream>
#include <assert.h>
#include <DMRG.hpp>
#include <ExactDiagonalization.hpp>

using namespace Eigen;
using namespace spinchain;


// with a bond dimension that can hold any state, DMRG must agree with exact diagonalization
int test_matches_exact(size_t chainlength, double spin, std::vector<double> couplings, std::vector<double> field)
{
        SpinChain model(chainlength, couplings, spin, field);

        ExactDiagonalization exact(model);
        double exactEnergy = exact.solve().energy;

        DMRGOptions options;
        options.chiMax = 256;
        options.truncCutoff = 1e-14;
        options.energyTol = 1e-11;
        options.lanczosTol = 1e-10;
        options.lanczosIterations = 20;
        options.maxSweeps = 40;
        DMRG dmrg(model, options);
        SolveResult result = dmrg.solve();

        assert( result.converged );
        assert( std::abs(result.energy - exactEnergy) < 1e-9 );
        assert( (dmrg.magnetization() - exact.magnetization()).norm() < 1e-6 );

        std::vector<std::complex<double>> dmrgSz = dmrg.localExpectations(model.site.Sz);
        std::vector<std::complex<double>> exactSz = exact.localExpectations(model.site.Sz);
        for (size_t i=0; i<chainlength; i++){
                assert( std::abs(dmrgSz[i] - exactSz[i]) < 1e-6 );
        }
        return 0;
}


// a singlet dimer has one bit of entanglement, which is ln 2
int test_entropy_of_dimer()
{
        DMRG dmrg(SpinChain(2, {1.0, 1.0, 1.0}, 0.5));
        assert( std::abs(dmrg.solve().energy + 0.75) < 1e-12 );
        assert( std::abs(dmrg.entanglementEntropy(0) - std::log(2.0)) < 1e-9 );
        return 0;
}


// a product state has no entanglement, and needs no truncation
int test_product_state()
{
        size_t chainlength = 12;
        DMRG dmrg(SpinChain(chainlength, {0.0, 0.0, 0.0}, 1.0, {0.0, 0.0, 2.0}));
        assert( std::abs(dmrg.solve().energy - (-2.0*chainlength)) < 1e-10 );
        for (size_t bond=0; bond+1<chainlength; bond++){
                assert( dmrg.entanglementEntropy(bond) < 1e-8 );
        }
        assert( dmrg.groundState().maxBondDim() <= 2 );
        return 0;
}


// results are only available after solve()
int test_requires_solve()
{
        DMRG dmrg(SpinChain(4, {1.0, 1.0, 1.0}, 0.5));
        bool threw = false;
        try {
                dmrg.magnetization();
        } catch (const std::logic_error &){
                threw = true;
        }
        assert( threw );
        return 0;
}


int main(){
        std::vector<double> xyz = {1.0, 0.7, 0.4};
        std::vector<double> tilted = {0.3, 0.2, 0.5};

        for (size_t chainlength : {2, 3, 4, 6, 8, 10}){
                test_matches_exact(chainlength, 0.5, xyz, tilted);
        }
        for (size_t chainlength : {2, 3, 4, 5, 6}){
                test_matches_exact(chainlength, 1.0, xyz, tilted);
        }
        test_matches_exact(4, 1.5, xyz, tilted);

        test_entropy_of_dimer();
        test_product_state();
        test_requires_solve();
}
