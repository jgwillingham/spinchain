import numpy as np
import pytest

import spinchain
from spinchain import DMRG, DMRGOptions, ExactDiagonalization, Spin, SpinChain

HEISENBERG = (1.0, 1.0, 1.0)


def solved_dmrg(length=8, **options):
    solver = DMRG(SpinChain(length, HEISENBERG, 0.5), DMRGOptions(chi_max=32, **options))
    return solver, solver.solve()


def test_spin_operators():
    site = Spin(0.5)
    assert site.dim == 2
    np.testing.assert_allclose(site.sz, np.diag([0.5, -0.5]))
    np.testing.assert_allclose(site.sx @ site.sy - site.sy @ site.sx, 1j * site.sz)


def test_model_attributes():
    model = SpinChain(4, [1.0, 2.0, 3.0], 1.0, field=(0.0, 0.0, 0.5))
    assert model.length == 4
    assert model.couplings == [1.0, 2.0, 3.0]
    assert model.field == [0.0, 0.0, 0.5]
    assert model.site.dim == 3


def test_exact_dimer_is_a_singlet():
    solver = ExactDiagonalization(SpinChain(2, HEISENBERG, 0.5))
    result = solver.solve()
    assert solver.dim == 4
    assert result.converged
    assert result.energy == pytest.approx(-0.75)
    np.testing.assert_allclose(solver.spectrum, [-0.75, 0.25, 0.25, 0.25], atol=1e-12)
    np.testing.assert_allclose(solver.magnetization(), np.zeros(3), atol=1e-12)


def test_exact_arrays_are_consistent():
    solver = ExactDiagonalization(SpinChain(4, (1.0, 0.5, 2.0), 0.5, (0.1, 0.0, 0.2)))
    solver.solve()
    H, states = solver.hamiltonian, solver.eigenstates
    np.testing.assert_allclose(H @ states, states @ np.diag(solver.spectrum), atol=1e-10)


def test_arrays_outlive_their_solver():
    solver = ExactDiagonalization(SpinChain(3, HEISENBERG, 0.5))
    solver.solve()
    spectrum = solver.spectrum
    del solver
    assert spectrum.shape == (8,)


def test_dmrg_matches_exact():
    model = SpinChain(8, (1.0, 0.8, 1.2), 0.5, (0.1, 0.0, 0.3))
    exact = ExactDiagonalization(model)
    dmrg = DMRG(model, DMRGOptions(chi_max=32))
    result = dmrg.solve()
    assert result.converged
    assert result.energy == pytest.approx(exact.solve().energy, abs=1e-7)
    np.testing.assert_allclose(dmrg.magnetization(), exact.magnetization(), atol=1e-5)


def test_dmrg_diagnostics():
    solver, _ = solved_dmrg()
    assert solver.max_truncation_error() < 1e-8
    assert solver.entanglement_entropy(3) > 0.0

    state = solver.ground_state()
    assert state.length == 8
    assert state.max_bond_dim <= 32
    assert state.entanglement_entropy(3) == solver.entanglement_entropy(3)
    assert state.tensor(0).shape[:2] == (1, 2)
    assert state.tensor(7).shape[1:] == (2, 1)


def test_site_tensors_contract_to_a_normalized_state():
    state = solved_dmrg(length=6)[0].ground_state()
    psi = state.tensor(0)
    for i in range(1, state.length):
        psi = np.tensordot(psi, state.tensor(i), axes=(-1, 0))
    assert np.vdot(psi, psi).real == pytest.approx(1.0)


def test_local_expectations_is_an_array():
    solver, _ = solved_dmrg()
    values = solver.local_expectations(Spin(0.5).sz)
    assert isinstance(values, np.ndarray)
    assert values.shape == (8,)


def test_options_are_keyword_only_with_library_defaults():
    options = DMRGOptions(chi_max=16)
    assert options.chi_max == 16
    assert options.max_sweeps == 30
    assert options.seed == 1
    options.seed = 7
    assert options.seed == 7
    with pytest.raises(TypeError):
        DMRGOptions(16)


@pytest.mark.parametrize("couplings", [(1.0, 1.0), (1.0, 1.0, 1.0, 1.0), ()])
def test_couplings_must_have_three_values(couplings):
    with pytest.raises(TypeError):
        SpinChain(4, couplings, 0.5)


def test_field_must_have_three_values():
    with pytest.raises(TypeError):
        SpinChain(4, HEISENBERG, 0.5, (1.0,))


def test_invalid_spin_raises():
    with pytest.raises(ValueError):
        Spin(0.3)
    with pytest.raises(ValueError):
        SpinChain(4, HEISENBERG, -1.0)


@pytest.mark.parametrize("make", [
    lambda: ExactDiagonalization(SpinChain(3, HEISENBERG, 0.5)),
    lambda: DMRG(SpinChain(3, HEISENBERG, 0.5)),
])
def test_observables_require_solve(make):
    solver = make()
    with pytest.raises(RuntimeError):
        solver.magnetization()
    with pytest.raises(RuntimeError):
        solver.local_expectations(Spin(0.5).sz)


@pytest.mark.parametrize("op", [np.eye(3), np.eye(2)[:1], np.zeros(2), np.zeros((0, 0))])
def test_wrong_operator_shape_raises(op):
    for solver in (ExactDiagonalization(SpinChain(3, HEISENBERG, 0.5)),
                   DMRG(SpinChain(3, HEISENBERG, 0.5))):
        solver.solve()
        with pytest.raises((ValueError, TypeError)):
            solver.local_expectations(op)


def test_dmrg_rejects_bad_arguments():
    with pytest.raises(ValueError):
        DMRG(SpinChain(1, HEISENBERG, 0.5))
    with pytest.raises(ValueError):
        DMRG(SpinChain(4, HEISENBERG, 0.5), DMRGOptions(chi_max=0))


@pytest.mark.parametrize("bond", [7, 8, 100])
def test_bond_out_of_range_raises(bond):
    solver, _ = solved_dmrg()
    with pytest.raises((IndexError, RuntimeError)):
        solver.entanglement_entropy(bond)


def test_site_out_of_range_raises():
    with pytest.raises(IndexError):
        solved_dmrg()[0].ground_state().tensor(8)


def test_overlong_exact_chain_raises():
    with pytest.raises(ValueError):
        ExactDiagonalization(SpinChain(40, HEISENBERG, 0.5))


def test_unaffordable_exact_chain_raises_memory_error():
    with pytest.raises(MemoryError):
        ExactDiagonalization(SpinChain(24, HEISENBERG, 0.5)).solve()


def test_non_finite_couplings_do_not_kill_the_interpreter():
    with pytest.raises(RuntimeError):
        DMRG(SpinChain(4, (float("nan"), 1.0, 1.0), 0.5)).solve()


def test_package_exports():
    assert set(spinchain.__all__) <= set(dir(spinchain))
