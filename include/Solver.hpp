#pragma once

#include <complex>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <SpinChain.hpp>

namespace spinchain{

using namespace Eigen;


/*! \brief Outcome of a ground state calculation
 */
struct SolveResult {
    double energy;      // ground state energy
    bool converged;     // false if the method stopped before meeting its tolerance
    size_t iterations;  // method-dependent: e.g. number of sweeps for DMRG
};


/*! \brief Interface shared by all methods for solving a SpinChain.
 *
 *  This covers only what every method can provide. Anything specific to one
 *  method, such as the full spectrum from exact diagonalization, belongs on the
 *  concrete class. To add a method, derive from Solver.
 */
class Solver {
    public:
        virtual ~Solver() = default;

        // compute the ground state
        virtual SolveResult solve() = 0;

        // returns the ground state expected value of op on each site; requires solve()
        virtual std::vector<std::complex<double>> localExpectations(const MatrixXcd& op) const = 0;

        // returns the expected value of the sum of spins as a 3D vector; requires solve()
        Vector3cd magnetization() const;

    protected:
        explicit Solver(const SpinChain& model);

        SpinChain model;
};

} // namespace spinchain
