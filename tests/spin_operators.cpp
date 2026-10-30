#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <eigen3/Eigen/Dense>
#include <Spin.hpp>

using namespace Eigen;
using namespace spinchain;
using namespace std::literals;



int test_spin_algebra(double spin)
{
        Spin s = Spin(spin);
        int dim = s.dim;
        assert( dim == int(2*spin + 1 + 0.5) );

        // [Sx, Sy] = i Sz
        MatrixXcd commutator = s.Sx*s.Sy - s.Sy*s.Sx;
        assert( commutator.isApprox( 1.0i*s.Sz, 1e-14 ) );

        // S^2 = S(S+1) * identity
        MatrixXcd S2 = s.Sx*s.Sx + s.Sy*s.Sy + s.Sz*s.Sz;
        assert( S2.isApprox( spin*(spin+1)*MatrixXcd::Identity(dim,dim), 1e-14 ) );
        return 0;
}



int test_invalid_spin()
{
        bool threw = false;
        try { Spin s = Spin(0.7); } catch (const std::invalid_argument&) { threw = true; }
        assert( threw );
        return 0;
}



int main(){
	test_spin_algebra(0.5);
	test_spin_algebra(1.0);
	test_spin_algebra(1.5);
	test_invalid_spin();
}
