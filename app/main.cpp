#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <SpinChain.hpp>


using namespace spinchain;

// the Hamiltonian is a dense complex matrix, so refuse sizes that would need gigabytes
static const double MAX_DIM = 8192;

static void printUsage(const char *prog, std::ostream &out)
{
	out << "usage: " << prog << " -L LENGTH [options]\n"
	    << "\n"
	    << "Exact diagonalization of the open spin-S XYZ chain\n"
	    << "  H = sum_i Jx Sx_i Sx_{i+1} + Jy Sy_i Sy_{i+1} + Jz Sz_i Sz_{i+1}\n"
	    << "\n"
	    << "options:\n"
	    << "  -L, --length N   number of sites (required, integer >= 2)\n"
	    << "  -s, --spin S     spin magnitude: 0.5, 1, 1.5, ...        (default 0.5)\n"
	    << "  -x, --jx J       Jx coupling                             (default 1)\n"
	    << "  -y, --jy J       Jy coupling                             (default 1)\n"
	    << "  -z, --jz J       Jz coupling                             (default 1)\n"
	    << "  -h, --help       show this message\n"
	    << "\n"
	    << "Values may be given as '--length 4' or '--length=4'.\n"
	    << "The Hilbert space dimension is (2S+1)^L, so keep chains short.\n"
	    << "\n"
	    << "examples:\n"
	    << "  " << prog << " -L 4                        spin-1/2 Heisenberg chain\n"
	    << "  " << prog << " -L 3 -s 1 --jz 0.5          anisotropic spin-1 chain\n"
	    << "  " << prog << " -L 6 -x 0 -y 0 -z 1         Ising chain\n";
}

// parses the whole string as a double, throwing if there is leftover text
static double parseNumber(const std::string &opt, const std::string &text)
{
	size_t used = 0;
	double value;
	try {
		value = std::stod(text, &used);
	} catch (const std::exception &) {
		throw std::invalid_argument("invalid value '" + text + "' for " + opt);
	}
	if (used != text.size() || !std::isfinite(value)){
		throw std::invalid_argument("invalid value '" + text + "' for " + opt);
	}
	return value;
}

int main(int argc, char **argv)
{
	size_t length = 0;
	bool haveLength = false;
	double spin = 0.5;
	std::vector<double> couplings = {1.0, 1.0, 1.0};

	try {
		for (int i = 1; i < argc; i++){
			std::string arg = argv[i];

			if (arg == "-h" || arg == "--help"){
				printUsage(argv[0], std::cout);
				return 0;
			}

			// split "--opt=value" into option and inline value
			std::string opt = arg, value;
			bool hasValue = false;
			size_t eq = arg.find('=');
			if (arg.rfind("--", 0) == 0 && eq != std::string::npos){
				opt = arg.substr(0, eq);
				value = arg.substr(eq + 1);
				hasValue = true;
			}

			if (opt != "-L" && opt != "--length" && opt != "-s" && opt != "--spin" &&
			    opt != "-x" && opt != "--jx" && opt != "-y" && opt != "--jy" &&
			    opt != "-z" && opt != "--jz"){
				throw std::invalid_argument("unknown argument '" + arg + "'");
			}

			// otherwise the value is the next argument
			if (!hasValue){
				if (i + 1 >= argc){
					throw std::invalid_argument("missing value for " + opt);
				}
				value = argv[++i];
			}

			if (opt == "-L" || opt == "--length"){
				double n = parseNumber(opt, value);
				if (n != std::floor(n) || n < 2){
					throw std::invalid_argument("chain length must be an integer >= 2");
				}
				length = static_cast<size_t>(n);
				haveLength = true;
			} else if (opt == "-s" || opt == "--spin"){
				spin = parseNumber(opt, value);
				double twoS = 2.0 * spin;
				if (spin <= 0 || std::abs(twoS - std::round(twoS)) > 1e-12){
					throw std::invalid_argument("spin must be a positive multiple of 1/2");
				}
			} else if (opt == "-x" || opt == "--jx"){
				couplings[0] = parseNumber(opt, value);
			} else if (opt == "-y" || opt == "--jy"){
				couplings[1] = parseNumber(opt, value);
			} else {
				couplings[2] = parseNumber(opt, value);
			}
		}

		if (!haveLength){
			throw std::invalid_argument("chain length is required (-L)");
		}

		double dim = std::pow(std::round(2.0 * spin) + 1.0, static_cast<double>(length));
		if (dim > MAX_DIM){
			throw std::invalid_argument("Hilbert space dimension " + std::to_string(dim) +
			    " exceeds the limit of " + std::to_string(static_cast<int>(MAX_DIM)) +
			    "; use a shorter chain or smaller spin");
		}
	} catch (const std::invalid_argument &e){
		std::cerr << argv[0] << ": " << e.what() << "\n\n";
		printUsage(argv[0], std::cerr);
		return 1;
	}

	std::cout << "Heisenberg XYZ model: L=" << length << ", S=" << spin << "\n";
	SpinChain heisenbergXYZ(length, couplings, spin);

	std::cout << "Jx = " << heisenbergXYZ.Jx << "\n";
	std::cout << "Jy = " << heisenbergXYZ.Jy << "\n";
	std::cout << "Jz = " << heisenbergXYZ.Jz << std::endl;

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

	return 0;
}
