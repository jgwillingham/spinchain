#include <iostream>
#include <vector>
#include <SpinHalfChain.hpp>

using namespace SpinChain;

int main(int argc, char **argv)
{
	std::vector<float> couplings;
	unsigned int length;

	length = 5;
	couplings = {1.0, 2.0, 3.0};

	SpinHalfChain heisenbergXYZ(length, couplings);

	for (int i=0; i< length; i++){
		std::cout << heisenbergXYZ.chain[i].Sz << std::endl;
	}

	

}
