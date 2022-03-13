#include <SpinHalf.hpp>

namespace SpinChain{

    using namespace std::literals;
    
    // define the spin-1/2 operators in the conventional z-eigenbasis
    // Each matrix is S_i = 0.5 * sigma_i
    
    const Matrix2cf SpinHalf::Sx  {{0.0f, 0.5f}, 
                                   {0.5f, 0.0f}};

    const Matrix2cf SpinHalf::Sy  {{0.0f,-0.5if}, 
                                   {0.5if, 0.0f}};

    const Matrix2cf SpinHalf::Sz   {{0.5f, 0.0f}, 
                                    {0.0f,-0.5f}};

}