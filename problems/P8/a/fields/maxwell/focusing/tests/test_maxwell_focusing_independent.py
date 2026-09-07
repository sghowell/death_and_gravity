"""Independent sampler integration and threshold reconstruction."""

import json
from pathlib import Path

from p8a_maxwell import verify as parent
from p8a_maxwell_focusing import controls, envelope, focusing, history, independent


def test_fraction_sampler_integrals_and_radical_margins_match_exact_API():
    report = independent.replay(json.loads(parent.REPORT.read_text()))
    assert report["cubic_moments"] == independent.serialize(history.moments())
    assert report["radical_margins"] == independent.serialize(envelope.radical_margins())


def test_fraction_threshold_and_complete_geometric_control_are_independent():
    report = independent.replay(json.loads(parent.REPORT.read_text()))
    constants, data = report["constants"], focusing.calibration()
    assert constants["total_cost"] == str(data["zero_source"]["envelope"]["total_cost"])
    assert constants["zero_source_margin"] == str(data["zero_source"]["strict_focusing_margin"])
    assert constants["sigma_one_margin"] == str(data["source_sigma_at_most_one"]["strict_focusing_margin"])
    assert report["geometric_countercontrol"] == independent.serialize(controls.calibration())


def test_fraction_engine_does_not_import_primary_formula_modules_or_sympy():
    source = Path(independent.__file__).read_text()
    assert "import sympy" not in source
    assert "from ." not in source
    assert "p8a_maxwell_focusing import" not in source
