#include <Environments.hpp>


namespace spinchain{


Tensor3 emptyEnvironment()
{
    Tensor3 environment(1,1,1);
    environment.setConstant(1.0);
    return environment;
}



Tensor3 growLeft(const Tensor3& left, const Tensor3& A, const Tensor4& W)
{
    Tensor3 bra = A.conjugate();

    // (bra, a, ket) * (ket, t, ketNew) = (bra, a, t, ketNew)
    Tensor4 withKet = left.contract(A, std::array<IndexPair,1>{IndexPair(2,0)});
    // ... * (a, b, s, t) = (bra, ketNew, b, s)
    Tensor4 withMpo = withKet.contract(W, std::array<IndexPair,2>{IndexPair(1,0), IndexPair(2,3)});
    // ... * (bra, s, braNew) = (ketNew, b, braNew)
    Tensor3 grown = withMpo.contract(bra, std::array<IndexPair,2>{IndexPair(0,0), IndexPair(3,1)});

    return grown.shuffle(std::array<Index,3>{2, 1, 0});
}



Tensor3 growRight(const Tensor3& right, const Tensor3& A, const Tensor4& W)
{
    Tensor3 bra = A.conjugate();

    // (ketNew, t, ket) * (bra, b, ket) = (ketNew, t, bra, b)
    Tensor4 withKet = A.contract(right, std::array<IndexPair,1>{IndexPair(2,2)});
    // ... * (a, b, s, t) = (ketNew, bra, a, s)
    Tensor4 withMpo = withKet.contract(W, std::array<IndexPair,2>{IndexPair(1,3), IndexPair(3,1)});
    // ... * (braNew, s, bra) = (ketNew, a, braNew)
    Tensor3 grown = withMpo.contract(bra, std::array<IndexPair,2>{IndexPair(1,2), IndexPair(3,1)});

    return grown.shuffle(std::array<Index,3>{2, 1, 0});
}



Tensor4 applyTwoSiteHamiltonian(const Tensor3& left, const Tensor4& W1, const Tensor4& W2,
                                const Tensor3& right, const Tensor4& theta)
{
    // (braL, a, ketL) * (ketL, s1, s2, ketR) = (braL, a, s1, s2, ketR)
    Eigen::Tensor<std::complex<double>,5> step1 =
        left.contract(theta, std::array<IndexPair,1>{IndexPair(2,0)});
    // ... * (a, b, s1Out, s1) = (braL, s2, ketR, b, s1Out)
    Eigen::Tensor<std::complex<double>,5> step2 =
        step1.contract(W1, std::array<IndexPair,2>{IndexPair(1,0), IndexPair(2,3)});
    // ... * (b, c, s2Out, s2) = (braL, ketR, s1Out, c, s2Out)
    Eigen::Tensor<std::complex<double>,5> step3 =
        step2.contract(W2, std::array<IndexPair,2>{IndexPair(3,0), IndexPair(1,3)});
    // ... * (braR, c, ketR) = (braL, s1Out, s2Out, braR)
    return step3.contract(right, std::array<IndexPair,2>{IndexPair(1,2), IndexPair(3,1)});
}



double energy(const MPS& state, const MPO& hamiltonian)
{
    Tensor3 environment = emptyEnvironment();
    for (size_t i=0; i<state.length(); i++){
        environment = growLeft(environment, state.tensor(i), hamiltonian.W[i]);
    }
    return environment(0,0,0).real();
}

} //namespace spinchain
