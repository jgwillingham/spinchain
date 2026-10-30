#pragma once

#include <complex>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <Tensors.hpp>

namespace spinchain{

using namespace Eigen;


/*! \brief Matrix product state for an open chain
 *
 *  The state is a product of site tensors with indices (bondLeft, physical, bondRight).
 *  It is always kept normalized and in mixed canonical form: tensors left of the
 *  orthogonality center are left isometries, and tensors right of it are right
 *  isometries. All the weight of the state then sits in the center tensor, so the
 *  rest of the chain acts as an orthonormal basis. This is what keeps DMRG stable.
 */
class MPS {
    public:
        /*! Random state with bond dimensions up to bondDim, and the center on site 0.
         *
         *  Being random, it overlaps every symmetry sector. This matters because
         *  DMRG cannot leave the sectors that its starting state occupies.
         */
        MPS(size_t length, Index physicalDim, Index bondDim, unsigned int seed);

        size_t length() const;
        size_t center() const;
        const Tensor3& tensor(size_t i) const;

        // the largest bond dimension in the chain
        Index maxBondDim() const;

        // moves the orthogonality center to the given site using QR decompositions
        void moveCenterTo(size_t site);

        // contracts sites i and i+1 into a tensor (bondLeft, physical1, physical2, bondRight)
        Tensor4 twoSiteTensor(size_t i) const;

        /*! Splits a two-site tensor back onto sites i and i+1, truncating the bond.
         *
         *  The center must be on site i or i+1, and it finishes on site i+1 if
         *  centerRight is true or on site i otherwise. Returns the discarded weight.
         */
        double setTwoSiteTensor(size_t i, const Tensor4& theta, Index maxBond, double cutoff,
                                bool centerRight);

        // returns the expected value of op on each site
        std::vector<std::complex<double>> localExpectations(const MatrixXcd& op) const;

        // returns the von Neumann entropy across the bond between sites i and i+1
        // from the most recent setTwoSiteTensor call on that bond
        double entanglementEntropy(size_t bond) const;

    private:
        void shiftCenterRight();
        void shiftCenterLeft();

        std::vector<Tensor3> tensors;
        std::vector<VectorXd> schmidtValues;    // singular values on each bond
        size_t centerSite;
};

} // namespace spinchain
