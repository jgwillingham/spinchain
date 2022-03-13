#include <SpinHalfChain.hpp>


namespace SpinChain{


/*! \brief Constructor for a chain of spin-1/2 objects
 *  
 *  When a SpinHalfChain object is declared, a vector of SpinHalf objects is assembled
 *  for the given chain length.
 */
SpinHalfChain::SpinHalfChain(unsigned int length, std::vector<float> couplings){
    // save couplings
    couplings = couplings;

    // assemble spin-1/2 chain
    for (int i=0; i < length; i++){
        SpinHalf chain[i];
    }
}


void SpinHalfChain::buildHamiltonian(){

    

}


}