import copy
import json

import pytest
from p8_composite_light import verify


def test_source_hashed_certificate_and_full_immutable_lineage():
    actual = verify.build_report()
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)
    assert actual["prior_S6_8_composite_cones_sha256"] == verify.PRIOR_SHA
    assert actual["claim"] == "P8-S6.9.COMPOSITE"


@pytest.mark.parametrize("key", ["domain", "retarded_and_local_representative",
                               "fixed_charge_zero_momentum", "strict_hierarchy_enclosures",
                               "not_established"])
def test_changed_domain_or_inflated_matching_scope_fails_replay(key):
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered[key] = "all ultraviolet parents and finite-band reductions excluded"
    with pytest.raises(ValueError, match="differs from exact replay"):
        verify.validate_report(altered, actual)


def test_all_state_and_time_dependent_omission_controls_replayed():
    assert verify.domain_controls()["zero_mass_inverse_representative"] == "rejected"
    actual = verify.build_report()
    assert actual["omission_controls"] == {
        "omitted_rotating_weight_heavy_forcing": "2",
        "omitted_homogeneous_heavy_data_light_forcing": "-2",
        "omitted_nonconstant_mass_derivatives_residual": "-12",
        "omitted_time_boundary_changes_momentum": "-1",
    }
