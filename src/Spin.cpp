#include <Spin.hpp>
#include <cmath>
#include <stdexcept>

namespace spinchain{

    using namespace std::literals;


/*! \brief Constructor for the spin-S operators
 *
 *  Sz is diagonal with entries m = S, ..., -S. Sx and Sy are built from the
 *  ladder operator S+, whose only nonzero elements are
 *  <m+1|S+|m> = sqrt(S(S+1) - m(m+1)).
 */
Spin::Spin(double spin) : spin(spin)
{
    // 2S must be a positive integer
    double twoS = 2.0*spin;
    if (spin <= 0.0 || std::abs(twoS - std::round(twoS)) > 1e-12){
        throw std::invalid_argument("spin must be a positive multiple of 1/2");
    }
    dim = static_cast<unsigned int>(std::round(twoS)) + 1;

    MatrixXcd Splus = MatrixXcd::Zero(dim,dim);
    Sz = MatrixXcd::Zero(dim,dim);

    for (unsigned int k=0; k<dim; k++){
        double m = spin - k;
        Sz(k,k) = m;
        if (k > 0){
            Splus(k-1,k) = std::sqrt(spin*(spin+1) - m*(m+1));
        }
    }
    MatrixXcd Sminus = Splus.adjoint();

    Sx = 0.5*(Splus + Sminus);
    Sy = -0.5i*(Splus - Sminus);
}

} //namespace spinchain
