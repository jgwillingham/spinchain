#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>
#include <Spin.hpp>

namespace spinchain{

using namespace Eigen;


/*! \brief One term of the two-site coupling: coupling * (left on site i) * (right on site i+1)
 */
struct BondTerm {
    double coupling;
    MatrixXcd left;
    MatrixXcd right;
};


/*! \brief Definition of a 1D spin-S model.
 *
 *  SpinChain describes an open spin-S chain with nearest-neighbor XYZ couplings
 *  in an external field. It only defines the model: solving it is the job of a
 *  Solver. Solvers build the Hamiltonian from bondTerms() and onsiteTerm(), so
 *  they all share one definition of the model.
 */
class SpinChain {
    public:
        /*! field is the external magnetic field (hx, hy, hz); it enters the
         *  Hamiltonian as the Zeeman term -sum_i (hx Sx_i + hy Sy_i + hz Sz_i).
         */
        SpinChain(size_t chainlength, std::vector<double> couplings, double spin,
                  std::vector<double> field = {0.0, 0.0, 0.0});

        size_t chainlength;             // Number of sites
        Spin site;                      // Single-site spin operators
        std::vector<double> couplings;  // Exchange couplings: Jx, Jy, Jz
        double Jx, Jy, Jz;
        std::vector<double> field;      // External magnetic field: hx, hy, hz
        double hx, hy, hz;

        // the two-site coupling, which is the same on every bond
        std::vector<BondTerm> bondTerms() const;

        // the single-site (Zeeman) term, which is the same on every site
        MatrixXcd onsiteTerm() const;
};


} // namespace spinchain
