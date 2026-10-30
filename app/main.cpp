#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <DMRG.hpp>
#include <ExactDiagonalization.hpp>


using namespace spinchain;

// exact diagonalization uses a dense complex matrix, so refuse sizes that would need gigabytes
static const double MAX_DIM = 8192;

static void printUsage(const char *prog, std::ostream &out)
{
	out << "usage: " << prog << " -L LENGTH [options]\n"
	    << "\n"
	    << "Ground state of the open spin-S XYZ chain\n"
	    << "  H = sum_i Jx Sx_i Sx_{i+1} + Jy Sy_i Sy_{i+1} + Jz Sz_i Sz_{i+1}\n"
	    << "      - sum_i (hx Sx_i + hy Sy_i + hz Sz_i)        (external field)\n"
	    << "\n"
	    << "options:\n"
	    << "  -L, --length N   number of sites (required, integer >= 2)\n"
	    << "  -s, --spin S     spin magnitude: 0.5, 1, 1.5, ...        (default 0.5)\n"
	    << "  -x, --jx J       Jx coupling                             (default 1)\n"
	    << "  -y, --jy J       Jy coupling                             (default 1)\n"
	    << "  -z, --jz J       Jz coupling                             (default 1)\n"
	    << "      --hx H       external field along x                  (default 0)\n"
	    << "      --hy H       external field along y                  (default 0)\n"
	    << "      --hz H       external field along z                  (default 0)\n"
	    << "  -m, --method M   solver: ed or dmrg                      (default ed)\n"
	    << "  -h, --help       show this message\n"
	    << "\n"
	    << "dmrg options:\n"
	    << "      --chi N      maximum bond dimension                  (default 100)\n"
	    << "      --sweeps N   maximum number of sweeps                (default 30)\n"
	    << "      --tol X      energy convergence tolerance            (default 1e-8)\n"
	    << "      --cutoff X   discarded weight at each truncation     (default 1e-12)\n"
	    << "\n"
	    << "Values may be given as '--length 4' or '--length=4'.\n"
	    << "\n"
	    << "ed finds the whole spectrum exactly, but the Hilbert space dimension is\n"
	    << "(2S+1)^L, so chains must be short. dmrg finds only the ground state, and\n"
	    << "reaches long chains. Its accuracy depends on --chi: compare a few values.\n"
	    << "\n"
	    << "examples:\n"
	    << "  " << prog << " -L 4                        spin-1/2 Heisenberg chain\n"
	    << "  " << prog << " -L 3 -s 1 --jz 0.5          anisotropic spin-1 chain\n"
	    << "  " << prog << " -L 6 -x 0 -y 0 -z 1         Ising chain\n"
	    << "  " << prog << " -L 6 -x 0 -y 0 --hx 0.5     transverse-field Ising chain\n"
	    << "  " << prog << " -L 100 -m dmrg --chi 64     long spin-1/2 Heisenberg chain\n";
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

// parses a number which must be a positive integer
static double parsePositiveInteger(const std::string &opt, const std::string &text)
{
	double n = parseNumber(opt, text);
	if (n != std::floor(n) || n < 1){
		throw std::invalid_argument(opt + " must be a positive integer");
	}
	return n;
}

static void printHeader(const std::string &method)
{
	std::cout << "\nSolving with " << method << " . . .\n" << std::endl;
}

static void runExactDiagonalization(const SpinChain &model)
{
	printHeader("exact diagonalization");
	ExactDiagonalization solver(model);
	solver.solve();

	std::cout << "The spectrum is :\n";
	std::cout << solver.spectrum << std::endl << std::endl;

	std::cout << "The ground state magnetization is :\n";
	std::cout << solver.magnetization() << std::endl;
}

static void runDMRG(const SpinChain &model, const DMRGOptions &options)
{
	printHeader("DMRG");
	DMRG solver(model, options);
	SolveResult result = solver.solve();

	std::cout << "\nThe ground state energy is : " << result.energy << "\n";
	std::cout << "Energy per site : " << result.energy/model.chainlength << "\n";
	std::cout << "Converged : " << (result.converged ? "yes" : "NO") << " after "
	          << result.iterations << " sweeps\n";
	std::cout << "Largest truncation error in the last sweep : "
	          << solver.maxTruncationError() << std::endl << std::endl;

	std::cout << "The ground state magnetization is :\n";
	std::cout << solver.magnetization() << std::endl << std::endl;

	std::cout << "The entanglement entropy of each bond is :\n";
	for (size_t bond = 0; bond + 1 < model.chainlength; bond++){
		std::cout << solver.entanglementEntropy(bond) << "\n";
	}
	std::cout << std::flush;
}

int main(int argc, char **argv)
{
	size_t length = 0;
	bool haveLength = false;
	double spin = 0.5;
	std::vector<double> couplings = {1.0, 1.0, 1.0};
	std::vector<double> field = {0.0, 0.0, 0.0};
	std::string method = "ed";
	DMRGOptions dmrgOptions;
	dmrgOptions.verbose = true;

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
			    opt != "-z" && opt != "--jz" && opt != "--hx" && opt != "--hy" &&
			    opt != "--hz" && opt != "-m" && opt != "--method" && opt != "--chi" &&
			    opt != "--sweeps" && opt != "--tol" && opt != "--cutoff"){
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
			} else if (opt == "-z" || opt == "--jz"){
				couplings[2] = parseNumber(opt, value);
			} else if (opt == "--hx"){
				field[0] = parseNumber(opt, value);
			} else if (opt == "--hy"){
				field[1] = parseNumber(opt, value);
			} else if (opt == "--hz"){
				field[2] = parseNumber(opt, value);
			} else if (opt == "-m" || opt == "--method"){
				method = value;
				if (method != "ed" && method != "dmrg"){
					throw std::invalid_argument("method must be ed or dmrg");
				}
			} else if (opt == "--chi"){
				dmrgOptions.chiMax = static_cast<Index>(parsePositiveInteger(opt, value));
			} else if (opt == "--sweeps"){
				dmrgOptions.maxSweeps = static_cast<size_t>(parsePositiveInteger(opt, value));
			} else if (opt == "--tol"){
				dmrgOptions.energyTol = parseNumber(opt, value);
			} else {
				dmrgOptions.truncCutoff = parseNumber(opt, value);
				if (dmrgOptions.truncCutoff < 0){
					throw std::invalid_argument("--cutoff must not be negative");
				}
			}
		}

		if (!haveLength){
			throw std::invalid_argument("chain length is required (-L)");
		}

		double dim = std::pow(std::round(2.0 * spin) + 1.0, static_cast<double>(length));
		if (method == "ed" && dim > MAX_DIM){
			throw std::invalid_argument("Hilbert space dimension " + std::to_string(dim) +
			    " exceeds the limit of " + std::to_string(static_cast<int>(MAX_DIM)) +
			    "; use a shorter chain or smaller spin, or --method dmrg");
		}
	} catch (const std::invalid_argument &e){
		std::cerr << argv[0] << ": " << e.what() << "\n\n";
		printUsage(argv[0], std::cerr);
		return 1;
	}

	std::cout << "Heisenberg XYZ model: L=" << length << ", S=" << spin << "\n";
	SpinChain heisenbergXYZ(length, couplings, spin, field);

	std::cout << "Jx = " << heisenbergXYZ.Jx << "\n";
	std::cout << "Jy = " << heisenbergXYZ.Jy << "\n";
	std::cout << "Jz = " << heisenbergXYZ.Jz << "\n";
	std::cout << "h  = (" << heisenbergXYZ.hx << ", " << heisenbergXYZ.hy << ", "
	          << heisenbergXYZ.hz << ")" << std::endl;

	if (method == "dmrg"){
		runDMRG(heisenbergXYZ, dmrgOptions);
	} else {
		runExactDiagonalization(heisenbergXYZ);
	}

	return 0;
}
