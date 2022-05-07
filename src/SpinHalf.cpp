#include <SpinHalf.hpp>

namespace SpinChain{

    using namespace std::literals;
    
    // define the spin-1/2 operators in the conventional z-eigenbasis
    // Each matrix is S_i = 0.5 * sigma_i
    
    const Matrix2cd SpinHalf::Sx  {{0.0, 0.5}, 
                                   {0.5, 0.0}};

    const Matrix2cd SpinHalf::Sy  {{0.0,-0.5i}, 
                                   {0.5i, 0.0}};

    const Matrix2cd SpinHalf::Sz   {{0.5, 0.0}, 
                                    {0.0,-0.5}};

}
