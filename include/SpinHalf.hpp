#include <eigen3/Eigen/Dense>

namespace SpinChain{

using namespace Eigen;


/*! \brief Class for spin-1/2 operators
 *
 *	SpinHalf simply builds a spin-1/2 operator in the standard z-eigenbasis out of
 *  the standard Pauli matrices. Natural units are used (hbar=1).
 */
class SpinHalf 
{
	public:
		SpinHalf(){};
		static const Matrix2cd Sx;
		static const Matrix2cd Sy;
		static const Matrix2cd Sz;	
};

} //namespace SpinChain
