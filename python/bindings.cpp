#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>
#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <DMRG.hpp>
#include <ExactDiagonalization.hpp>


namespace nb = nanobind;
using namespace nb::literals;
using namespace spinchain;

using Couplings = std::array<double,3>;
using ComplexTensor = nb::ndarray<nb::numpy, std::complex<double>, nb::f_contig>;


/*! \brief A solver whose arguments are checked before they reach the library.
 *
 *  The library keeps assert() enabled, and a failed Eigen assert aborts the process.
 *  From Python that would kill the interpreter, so anything which would trip one has
 *  to be raised as an exception here instead.
 */
template <typename Method>
struct Checked : Method {
    using Method::Method;

    VectorXcd expectations(const MatrixXcd& op) const
    {
        Index dim = this->model.site.dim;
        if (op.rows() != dim || op.cols() != dim){
            throw std::invalid_argument("op must be a " + std::to_string(dim) + "x" +
                                        std::to_string(dim) + " matrix");
        }
        std::vector<std::complex<double>> values = this->localExpectations(op);
        return Map<const VectorXcd>(values.data(), static_cast<Index>(values.size()));
    }
};


// ExactDiagonalization stores its dimension as an unsigned int
static void checkExactDimension(const SpinChain& model)
{
    double dim = std::pow(static_cast<double>(model.site.dim), static_cast<double>(model.chainlength));
    if (dim > std::numeric_limits<unsigned int>::max()){
        throw std::invalid_argument("chain is too long for exact diagonalization; use DMRG");
    }
}


static ComplexTensor siteTensor(const MPS& state, size_t i)
{
    if (i >= state.length()){
        throw std::out_of_range("tensor: site is outside the chain");
    }

    // the array owns a copy, since the next solve replaces the tensors in the state
    Tensor3* copy = new Tensor3(state.tensor(i));
    nb::capsule owner(copy, [](void* p) noexcept { delete static_cast<Tensor3*>(p); });
    return ComplexTensor(copy->data(),
                         {static_cast<size_t>(copy->dimension(0)),
                          static_cast<size_t>(copy->dimension(1)),
                          static_cast<size_t>(copy->dimension(2))},
                         owner);
}


// the interface shared by all solvers; Eigen values are returned by copy for the same reason
template <typename Method>
static void bindSolver(nb::class_<Method>& solver)
{
    solver.def("solve", [](Method& self){ return self.solve(); },
               nb::call_guard<nb::gil_scoped_release>(),
               "Compute the ground state, releasing the GIL while it runs.")
          .def("local_expectations", &Method::expectations, "op"_a,
               "Ground state expected value of op on each site.")
          .def("magnetization", [](const Method& self){ return self.magnetization(); },
               "Ground state expected value of the sum of spins as a 3-vector.");
}


NB_MODULE(_core, m)
{
    m.doc() = "Quantum spin chains using exact diagonalization and DMRG";

    nb::class_<Spin>(m, "Spin", "Spin-S operators in the z-eigenbasis, ordered m = S, ..., -S.")
        .def(nb::init<double>(), "spin"_a)
        .def_ro("spin", &Spin::spin)
        .def_ro("dim", &Spin::dim)
        .def_prop_ro("sx", [](const Spin& self) -> MatrixXcd { return self.Sx; })
        .def_prop_ro("sy", [](const Spin& self) -> MatrixXcd { return self.Sy; })
        .def_prop_ro("sz", [](const Spin& self) -> MatrixXcd { return self.Sz; });

    nb::class_<SpinChain>(m, "SpinChain", "An open spin-S chain with XYZ couplings in a field.")
        .def("__init__",
             [](SpinChain* self, size_t length, Couplings couplings, double spin, Couplings field){
                 new (self) SpinChain(length, {couplings.begin(), couplings.end()}, spin,
                                      {field.begin(), field.end()});
             },
             "length"_a, "couplings"_a, "spin"_a, "field"_a = Couplings{0.0, 0.0, 0.0})
        .def_ro("length", &SpinChain::chainlength)
        .def_ro("site", &SpinChain::site)
        .def_ro("couplings", &SpinChain::couplings)
        .def_ro("field", &SpinChain::field);

    nb::class_<SolveResult>(m, "SolveResult", "Outcome of a ground state calculation.")
        .def_ro("energy", &SolveResult::energy)
        .def_ro("converged", &SolveResult::converged)
        .def_ro("iterations", &SolveResult::iterations)
        .def("__repr__", [](const SolveResult& self){
            return "SolveResult(energy=" + std::to_string(self.energy) +
                   ", converged=" + (self.converged ? "True" : "False") +
                   ", iterations=" + std::to_string(self.iterations) + ")";
        });

    const DMRGOptions defaults;
    nb::class_<DMRGOptions>(m, "DMRGOptions", "Settings for a DMRG calculation.")
        .def("__init__",
             [](DMRGOptions* self, Index chiMax, size_t maxSweeps, double energyTol,
                double truncCutoff, size_t lanczosIterations, double lanczosTol,
                unsigned int seed, bool verbose){
                 new (self) DMRGOptions{chiMax, maxSweeps, energyTol, truncCutoff,
                                        lanczosIterations, lanczosTol, seed, verbose};
             },
             nb::kw_only(),
             "chi_max"_a = defaults.chiMax,
             "max_sweeps"_a = defaults.maxSweeps,
             "energy_tol"_a = defaults.energyTol,
             "trunc_cutoff"_a = defaults.truncCutoff,
             "lanczos_iterations"_a = defaults.lanczosIterations,
             "lanczos_tol"_a = defaults.lanczosTol,
             "seed"_a = defaults.seed,
             "verbose"_a = defaults.verbose)
        .def_rw("chi_max", &DMRGOptions::chiMax)
        .def_rw("max_sweeps", &DMRGOptions::maxSweeps)
        .def_rw("energy_tol", &DMRGOptions::energyTol)
        .def_rw("trunc_cutoff", &DMRGOptions::truncCutoff)
        .def_rw("lanczos_iterations", &DMRGOptions::lanczosIterations)
        .def_rw("lanczos_tol", &DMRGOptions::lanczosTol)
        .def_rw("seed", &DMRGOptions::seed)
        .def_rw("verbose", &DMRGOptions::verbose);

    nb::class_<MPS>(m, "MPS", "Matrix product state for an open chain.")
        .def_prop_ro("length", &MPS::length)
        .def_prop_ro("center", &MPS::center)
        .def_prop_ro("max_bond_dim", &MPS::maxBondDim)
        .def("tensor", &siteTensor, "site"_a,
             "Copy of a site tensor with indices (bond_left, physical, bond_right).")
        .def("entanglement_entropy", &MPS::entanglementEntropy, "bond"_a);

    using PyExact = Checked<ExactDiagonalization>;
    nb::class_<PyExact> exact(m, "ExactDiagonalization", "Full spectrum of a short chain.");
    exact.def("__init__",
              [](PyExact* self, const SpinChain& model){
                  checkExactDimension(model);
                  new (self) PyExact(model);
              },
              "model"_a)
         .def_ro("dim", &PyExact::dim)
         .def("build_hamiltonian", &PyExact::buildHamiltonian)
         .def_prop_ro("hamiltonian", [](const PyExact& self) -> MatrixXcd { return self.hamiltonian; })
         .def_prop_ro("spectrum", [](const PyExact& self) -> VectorXd { return self.spectrum.real(); })
         .def_prop_ro("eigenstates", [](const PyExact& self) -> MatrixXcd { return self.eigenstates; });
    bindSolver(exact);

    using PyDMRG = Checked<DMRG>;
    nb::class_<PyDMRG> dmrg(m, "DMRG", "Two-site finite DMRG for the ground state.");
    dmrg.def(nb::init<const SpinChain&, DMRGOptions>(), "model"_a, "options"_a = DMRGOptions())
        .def("entanglement_entropy", &PyDMRG::entanglementEntropy, "bond"_a,
             "Von Neumann entropy across the bond between sites i and i+1.")
        .def("max_truncation_error", &PyDMRG::maxTruncationError,
             "Largest discarded weight of any truncation in the last sweep.")
        .def("ground_state", &PyDMRG::groundState, nb::rv_policy::reference_internal);
    bindSolver(dmrg);
}
