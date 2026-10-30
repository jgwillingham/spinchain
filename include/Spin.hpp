#pragma once

#include <eigen3/Eigen/Dense>

namespace spinchain{

using namespace Eigen;


/*! \brief Class for spin-S operators
 *
 *	Spin builds the spin operators Sx, Sy, Sz for a single spin of magnitude S
 *  (S = 1/2, 1, 3/2, ...) in the standard z-eigenbasis, ordered m = S, S-1, ..., -S.
 *  Natural units are used (hbar=1).
 */
class Spin
{
	public:
		Spin(double spin);

		double spin;		// spin magnitude S
		unsigned int dim;	// local Hilbert space dimension 2S+1
		MatrixXcd Sx;
		MatrixXcd Sy;
		MatrixXcd Sz;
};

} //namespace spinchain
