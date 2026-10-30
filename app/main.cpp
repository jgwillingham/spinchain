#include <iostream>
#include <vector>
#include <SpinChain.hpp>


using namespace spinchain;

int main(int argc, char **argv)
{
	if (argc < 2){
		std::cerr << "usage: " << argv[0] << " chainlength [spin]\n";
		return 1;
	}

	std::vector<double> couplings;
	unsigned int length;
	double spin;

	length = atoi(argv[1]);
	spin = (argc > 2) ? atof(argv[2]) : 0.5;
	couplings = {0.0, 0.0, 1.0};

	std::cout << "Heisenbery XYZ model:\n";
	SpinChain heisenbergXYZ(length, couplings, spin);

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
	std::cout << "The ground state magnetization is :\n";
        Vector3cd mag = heisenbergXYZ.getMagnetization(groundstate);
        std::cout << mag << std::endl;


}
