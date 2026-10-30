#include <MPS.hpp>
#include <SchmidtSplit.hpp>
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>


namespace spinchain{


static void normalize(Tensor3& tensor)
{
    Map<VectorXcd>(tensor.data(), tensor.size()).normalize();
}


/*
 * Extends the overlap <psi|op|psi> of everything left of a site across that site:
 *
 *      .---- conj(A) ---- bra
 *      |       |
 *    block     op (if given)
 *      |       |
 *      .------ A -------- ket
 */
static Tensor2 extendLeft(const Tensor2& block, const Tensor3& A, const Tensor2* op)
{
    // (bra, ket) * (ket, s, ketNew) = (bra, s, ketNew)
    Tensor3 partial = block.contract(A, std::array<IndexPair,1>{IndexPair(1,0)});
    if (op != nullptr){
        // (bra, s, ketNew) * (sOut, s) = (bra, ketNew, sOut), and then back to (bra, sOut, ketNew)
        Tensor3 acted = partial.contract(*op, std::array<IndexPair,1>{IndexPair(1,1)});
        partial = acted.shuffle(std::array<Index,3>{0, 2, 1});
    }
    // (braNew, ketNew)
    Tensor3 bra = A.conjugate();
    Tensor2 extended = bra.contract(partial, std::array<IndexPair,2>{IndexPair(0,0), IndexPair(1,1)});
    return extended;
}


// mirror image of extendLeft for everything right of a site, without an operator
static Tensor2 extendRight(const Tensor2& block, const Tensor3& A)
{
    // (ketNew, s, ket) * (bra, ket) = (ketNew, s, bra)
    Tensor3 partial = A.contract(block, std::array<IndexPair,1>{IndexPair(2,1)});
    // (braNew, ketNew)
    Tensor3 bra = A.conjugate();
    Tensor2 extended = bra.contract(partial, std::array<IndexPair,2>{IndexPair(2,2), IndexPair(1,1)});
    return extended;
}



MPS::MPS(size_t length, Index physicalDim, Index bondDim, unsigned int seed)
    : schmidtValues(length > 0 ? length - 1 : 0), centerSite(0)
{
    if (length == 0 || physicalDim < 1 || bondDim < 1){
        throw std::invalid_argument("MPS requires positive length and dimensions");
    }

    // a bond cannot usefully be larger than the Hilbert space on either side of it
    std::vector<Index> bonds(length+1, 1);
    for (size_t i=1; i<length; i++){
        bonds[i] = std::min(bondDim, bonds[i-1]*physicalDim);
    }
    for (size_t i=length-1; i>0; i--){
        bonds[i] = std::min(bonds[i], bonds[i+1]*physicalDim);
    }

    std::mt19937 generator(seed);
    std::normal_distribution<double> gaussian(0.0, 1.0);
    for (size_t i=0; i<length; i++){
        Tensor3 site(bonds[i], physicalDim, bonds[i+1]);
        for (Index k=0; k<site.size(); k++){
            site.data()[k] = std::complex<double>(gaussian(generator), gaussian(generator));
        }
        normalize(site);
        tensors.push_back(site);
    }

    // sweeping from the right turns every tensor it passes into a right isometry
    centerSite = length - 1;
    moveCenterTo(0);
}



size_t MPS::length() const
{
    return tensors.size();
}



size_t MPS::center() const
{
    return centerSite;
}



const Tensor3& MPS::tensor(size_t i) const
{
    return tensors.at(i);
}



Index MPS::maxBondDim() const
{
    Index largest = 1;
    for (const Tensor3& site : tensors){
        largest = std::max(largest, site.dimension(2));
    }
    return largest;
}



void MPS::shiftCenterRight()
{
    Tensor3& site = tensors[centerSite];
    Tensor3& next = tensors[centerSite+1];
    Index bondLeft = site.dimension(0);
    Index d = site.dimension(1);
    Index bondRight = site.dimension(2);

    // site = Q*R, where Q stays here as a left isometry and R moves to the next site
    HouseholderQR<MatrixXcd> qr(Map<const MatrixXcd>(site.data(), bondLeft*d, bondRight));
    Index k = std::min(bondLeft*d, bondRight);
    MatrixXcd Q = qr.householderQ()*MatrixXcd::Identity(bondLeft*d, k);
    MatrixXcd R = qr.matrixQR().topRows(k).triangularView<Upper>();

    Tensor3 moved = TensorMap<Tensor2>(R.data(), k, bondRight)
                        .contract(next, std::array<IndexPair,1>{IndexPair(1,0)});
    next = moved;
    site = TensorMap<Tensor3>(Q.data(), bondLeft, d, k);
    normalize(next);
    centerSite++;
}



void MPS::shiftCenterLeft()
{
    Tensor3& site = tensors[centerSite];
    Tensor3& previous = tensors[centerSite-1];
    Index bondLeft = site.dimension(0);
    Index d = site.dimension(1);
    Index bondRight = site.dimension(2);

    // site^dagger = Q*R, so site = R^dagger * Q^dagger, where Q^dagger stays here as a
    // right isometry and R^dagger moves to the previous site
    MatrixXcd adjoint = Map<const MatrixXcd>(site.data(), bondLeft, d*bondRight).adjoint();
    HouseholderQR<MatrixXcd> qr(adjoint);
    Index k = std::min(bondLeft, d*bondRight);
    MatrixXcd Q = qr.householderQ()*MatrixXcd::Identity(d*bondRight, k);
    MatrixXcd R = qr.matrixQR().topRows(k).triangularView<Upper>();
    MatrixXcd Rdagger = R.adjoint();
    MatrixXcd Qdagger = Q.adjoint();

    Tensor3 moved = previous.contract(TensorMap<Tensor2>(Rdagger.data(), bondLeft, k),
                                      std::array<IndexPair,1>{IndexPair(2,0)});
    previous = moved;
    site = TensorMap<Tensor3>(Qdagger.data(), k, d, bondRight);
    normalize(previous);
    centerSite--;
}



void MPS::moveCenterTo(size_t site)
{
    if (site >= tensors.size()){
        throw std::out_of_range("moveCenterTo: site is outside the chain");
    }
    while (centerSite < site){
        shiftCenterRight();
    }
    while (centerSite > site){
        shiftCenterLeft();
    }
}



Tensor4 MPS::twoSiteTensor(size_t i) const
{
    return tensors.at(i).contract(tensors.at(i+1), std::array<IndexPair,1>{IndexPair(2,0)});
}



double MPS::setTwoSiteTensor(size_t i, const Tensor4& theta, Index maxBond, double cutoff,
                             bool centerRight)
{
    if (centerSite != i && centerSite != i+1){
        throw std::logic_error("setTwoSiteTensor: the center must be on one of the two sites");
    }
    Index bondLeft = theta.dimension(0);
    Index d1 = theta.dimension(1);
    Index d2 = theta.dimension(2);
    Index bondRight = theta.dimension(3);

    SchmidtSplit split = schmidtSplit(Map<const MatrixXcd>(theta.data(), bondLeft*d1, d2*bondRight),
                                      maxBond, cutoff, centerRight);
    Index k = split.schmidtValues.size();

    tensors[i] = TensorMap<Tensor3>(split.left.data(), bondLeft, d1, k);
    tensors[i+1] = TensorMap<Tensor3>(split.right.data(), k, d2, bondRight);
    schmidtValues[i] = split.schmidtValues;
    centerSite = centerRight ? i+1 : i;

    return split.discardedWeight;
}



std::vector<std::complex<double>> MPS::localExpectations(const MatrixXcd& op) const
{
    size_t length = tensors.size();
    Tensor2 opTensor = TensorMap<const Tensor2>(op.data(), op.rows(), op.cols());

    // rightBlocks[i] is the overlap of everything right of site i
    std::vector<Tensor2> rightBlocks(length);
    rightBlocks[length-1] = Tensor2(1,1);
    rightBlocks[length-1].setConstant(1.0);
    for (size_t i=length-1; i>0; i--){
        rightBlocks[i-1] = extendRight(rightBlocks[i], tensors[i]);
    }

    std::vector<std::complex<double>> values;
    Tensor2 leftBlock(1,1);
    leftBlock.setConstant(1.0);
    for (size_t i=0; i<length; i++){
        Tensor2 withOp = extendLeft(leftBlock, tensors[i], &opTensor);
        Eigen::Tensor<std::complex<double>,0> value =
            withOp.contract(rightBlocks[i], std::array<IndexPair,2>{IndexPair(0,0), IndexPair(1,1)});
        values.push_back(value());

        Tensor2 extended = extendLeft(leftBlock, tensors[i], nullptr);
        leftBlock = extended;
    }
    return values;
}



double MPS::entanglementEntropy(size_t bond) const
{
    const VectorXd& s = schmidtValues.at(bond);
    if (s.size() == 0){
        throw std::logic_error("entanglementEntropy: this bond has not been optimized");
    }

    double entropy = 0.0;
    for (Index k=0; k<s.size(); k++){
        double weight = s(k)*s(k);
        if (weight > 0.0){
            entropy -= weight*std::log(weight);
        }
    }
    return entropy;
}

} //namespace spinchain
