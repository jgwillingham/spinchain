#include <SpinChain.hpp>


namespace spinchain{


/*! \brief Constructor for a chain of spin-S objects
 *
 *  The spin magnitude is shared by every site, so a single Spin object holds
 *  the local operators for the whole chain.
 */
SpinChain::SpinChain(size_t chainlength, std::vector<double> couplings, double spin,
                     std::vector<double> field)
    : chainlength(chainlength), site(spin)
{
    // save couplings
    this->couplings = couplings;
    Jx = couplings[0];
    Jy = couplings[1];
    Jz = couplings[2];

    // save external field
    this->field = field;
    hx = field[0];
    hy = field[1];
    hz = field[2];
}



std::vector<BondTerm> SpinChain::bondTerms() const
{
    return {{Jx, site.Sx, site.Sx},
            {Jy, site.Sy, site.Sy},
            {Jz, site.Sz, site.Sz}};
}



MatrixXcd SpinChain::onsiteTerm() const
{
    return -(hx*site.Sx + hy*site.Sy + hz*site.Sz);
}

} //namespace spinchain
