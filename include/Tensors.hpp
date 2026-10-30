#pragma once

#include <array>
#include <complex>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

namespace spinchain{

// Tensors store their first index contiguously (column major), so a tensor can be
// viewed as a matrix without copying by splitting its indices into two groups.
using Tensor2 = Eigen::Tensor<std::complex<double>, 2>;
using Tensor3 = Eigen::Tensor<std::complex<double>, 3>;
using Tensor4 = Eigen::Tensor<std::complex<double>, 4>;

// a pair of indices to sum over in Tensor::contract
using IndexPair = Eigen::IndexPair<Eigen::Index>;

} // namespace spinchain
