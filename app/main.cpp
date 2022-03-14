#include <iostream>
#include <vector>
#include <SpinHalfChain.hpp>
#include <eigen3/unsupported/Eigen/KroneckerProduct>


using namespace SpinChain;

int main(int argc, char **argv)
{
	std::vector<float> couplings;
	unsigned int length;

	length = atoi(argv[1]);
	couplings = {1.0, 2.0, 3.0};

	std::cout << "Heisenbery XYZ model:\n\n";

	SpinHalfChain heisenbergXYZ(length, couplings);

	heisenbergXYZ.buildHamiltonian();

	std::cout << heisenbergXYZ.hamiltonian << std::endl << std::endl;

	heisenbergXYZ.diagonalize();
	std::cout << "Diagonalizing . . .\n\n" << std::endl;
	std::cout << "The ground state energy is :" << std::endl;
	std::cout << heisenbergXYZ.spectrum[0] << std::endl;

	std::cout << "The ground state is :" << std::endl;
	std::cout << heisenbergXYZ.eigenstates.col(0) << std::endl;
	

}
