import copy
import json

import pytest
from p8_scattering import verify


def test_scattering_report_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_changed_majorant_and_uv_verdict_rejected():
    actual = verify.build_report()
    bad = copy.deepcopy(actual)
    bad['uniform_finite_time_tree_majorant']['sufficient_M_tau'] = '1'
    with pytest.raises(ValueError, match='differs'):
        verify.validate_report(bad, actual)
    bad = copy.deepcopy(actual)
    bad['uv_applicability']['positivity_verdict'] = 'PASSED'
    with pytest.raises(ValueError, match='differs'):
        verify.validate_report(bad, actual)
