#include <iostream>
#include <vector>
#include <SpinHalfChain.hpp>


using namespace SpinChain;

int main(int argc, char **argv)
{
	std::vector<double> couplings;
	unsigned int length;

	length = atoi(argv[1]);
	couplings = {1.0, 0.0, 0.0};

	std::cout << "Heisenbery XYZ model:\n";
	SpinHalfChain heisenbergXYZ(length, couplings);
	
	std::cout << heisenbergXYZ.Jx << std::endl;
	std::cout << heisenbergXYZ.Jy << std::endl;
	std::cout << heisenbergXYZ.Jz << std::endl;

	heisenbergXYZ.buildHamiltonian();
	std::cout << "Hamiltonian made\n";

	heisenbergXYZ.diagonalize();
	std::cout << "Diagonalizing . . .\n" << std::endl;

 	std::cout << "The spectrum is :\n";
	std::cout << heisenbergXYZ.spectrum << std::endl << std::endl;

	VectorXcd groundstate = heisenbergXYZ.eigenstates.col(0);

        
	//std::cout << "The magnetization is :\n";
        //Vector3cf mag = heisenbergXYZ.getMagnetization(groundstate);
        //std::cout << mag << std::endl;
	

}
