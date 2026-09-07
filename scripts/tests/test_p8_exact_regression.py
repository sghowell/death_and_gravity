"""Independent boundary and restoration controls for the opt-in P8 runner."""

import importlib.util
from collections import Counter
from pathlib import Path

import pytest
from sympy.polys.domains import QQ, QQ_I, ZZ_I
from sympy.polys.rings import PolyElement, ring

SPEC = importlib.util.spec_from_file_location(
    "p8_exact_regression", Path(__file__).resolve().parents[1]/"p8_exact_regression.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


def test_reference_gaussian_tuples_for_all_unit_phases_and_rational_contents():
    assert runner.self_check() == {"original_tuple_comparisons": 128, "exact_descent": 128}


@pytest.mark.parametrize("domain", [QQ_I, ZZ_I])
def test_mixed_complex_coefficients_use_unchanged_algorithm(domain):
    _, x = ring("x", domain)
    imaginary = domain.dtype(0, 1)
    left, right = (x+imaginary)*(x+1), (x+imaginary)*(x+2)
    counts = Counter()
    assert runner.descended_gcd(left, right, counts) == runner.ORIGINAL_GCD(left, right)
    assert counts == {"mixed_fallback": 1}


@pytest.mark.parametrize("case", ["multivariate", "real_domain", "zero"])
def test_unaudited_domains_and_zero_operands_use_unchanged_algorithm(case):
    if case == "multivariate":
        _, x, y = ring("x,y", QQ_I)
        left, right = (x+y)*(x+1), (x+y)*(y+1)
    else:
        _, x = ring("x", QQ if case == "real_domain" else QQ_I)
        left, right = (x+1)*(x+2), (x+1)*(x+3)
        if case == "zero":
            right = right.ring.zero
    counts = Counter()
    assert runner.descended_gcd(left, right, counts) == runner.ORIGINAL_GCD(left, right)
    assert counts == {"domain_fallback": 1}


def test_context_restores_original_method_after_exception():
    original = PolyElement._gcd
    with pytest.raises(ValueError, match="intentional"), runner.exact_runner():
        assert PolyElement._gcd is not original
        raise ValueError("intentional")
    assert PolyElement._gcd is original


def test_other_version_is_not_silently_admitted(monkeypatch):
    monkeypatch.setattr(runner.sp, "__version__", "unreviewed")
    with pytest.raises(RuntimeError, match="audited only"):
        runner.self_check()
    with pytest.raises(RuntimeError, match="Unexpected"), runner.exact_runner():
        pass


def test_existing_replacement_is_not_overwritten(monkeypatch):
    monkeypatch.setattr(PolyElement, "_gcd", lambda left, right: None)
    with pytest.raises(RuntimeError, match="existing"), runner.exact_runner():
        pass


def test_bad_real_cofactor_result_fails_closed(monkeypatch):
    _, x = ring("x", QQ_I)
    monkeypatch.setattr(runner, "ORIGINAL_GCD", lambda left, right: (left, left, right))
    with pytest.raises(ArithmeticError, match="cofactor"):
        runner.descended_gcd((x+1)*(x+2), (x+1)*(x+3), Counter())


@pytest.mark.parametrize("arguments", [[], ["problems/P8", "-q"], ("problems/P8", "--collect-only")])
def test_default_collection_isolates_duplicate_frozen_basenames(arguments):
    before = list(arguments)
    assert runner.pytest_arguments(arguments) == ["--import-mode=importlib", *before]
    assert list(arguments) == before


@pytest.mark.parametrize("mode", ["prepend", "append", "importlib"])
@pytest.mark.parametrize("separate", [False, True])
def test_explicit_import_mode_is_preserved(mode, separate):
    option = ["--import-mode", mode] if separate else ["--import-mode="+mode]
    arguments = ["problems/P8", *option, "-q"]
    assert runner.pytest_arguments(arguments) == arguments
    assert runner.pytest_arguments(arguments) is not arguments


@pytest.mark.parametrize("exit_code", [0, 2])
def test_main_passes_isolated_arguments_and_restores_exact_runner(monkeypatch, exit_code):
    received = []
    original = PolyElement._gcd

    def fake_pytest(arguments):
        assert PolyElement._gcd is not original
        received.append(arguments)
        return exit_code

    monkeypatch.setattr(pytest, "main", fake_pytest)
    assert runner.main(["problems/P8", "-q"]) == exit_code
    assert received == [["--import-mode=importlib", "problems/P8", "-q"]]
    assert PolyElement._gcd is original
