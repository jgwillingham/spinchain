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
		static const Matrix2cf Sx;
		static const Matrix2cf Sy;
		static const Matrix2cf Sz;	
};

} //namespace SpinChain
