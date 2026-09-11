"""Disjoint four-point families and their already assigned counterterms."""

import sympy as s


def rows():
    return [
        {
            "id": "scalar_UV_finite",
            "refinements": 88,
            "owner": "S6.121",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_OS_insertion",
            "refinements": 64,
            "owner": "S6.153",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_double_bubble",
            "refinements": 24,
            "owner": "S6.152",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_wineglass",
            "refinements": 16,
            "owner": "S6.151",
            "status": "BOUNDED",
        },
        {"id": "scalar_Phi4_W0_F4", "owner": "S6.142", "status": "BOUNDED"},
        {"id": "gauge_Phi4", "owner": "S6.142", "status": "BOUNDED"},
        {"id": "scalar_Phi4_W1_F2", "owner": "S6.144", "status": "BOUNDED"},
        {"id": "scalar_Phi4_W2_F0", "owner": "S6.143", "status": "BOUNDED"},
    ]


def counterterms():
    return [
        {
            "id": "scalar_proper_interaction_forests",
            "owner": "scalar_192",
            "extra": False,
        },
        {
            "id": "scalar_inner_OS_mass_and_kinetic",
            "owner": "scalar_64",
            "extra": False,
        },
        {
            "id": "fermion_mass_kinetic_and_Yukawa_in_direct_rows",
            "owner": "S6.142",
            "extra": False,
        },
        {
            "id": "fermion_box_proper_reference_and_outer_MS",
            "owner": "S6.144",
            "extra": False,
        },
        {"id": "fermion_OS_insertion_and_outer_MS", "owner": "S6.143", "extra": False},
        {
            "id": "isolated_sigma_full_functional_conversion",
            "owner": "S6.150",
            "extra": False,
        },
        {
            "id": "full_regulated_first_and_second_Phi_map",
            "owner": "S6.156",
            "extra": True,
        },
        {
            "id": "first_coordinate_variation_of_finite_one_loop",
            "owner": "this_checkpoint",
            "extra": True,
        },
        {"id": "overall_local_Phi4_reference", "owner": "zero_b2", "extra": False},
        {
            "id": "Phi2_mass_source_and_external_residue",
            "owner": "S6.156_OS",
            "extra": False,
        },
    ]


def validate_rows(value, ct):
    if value != rows() or ct != counterterms():
        raise ValueError("The complete amplitude ownership ledger changed")
    return True


def data():
    J, G, Phi, M = s.symbols("J G Phi M", nonzero=True)
    eliminated = -((J + G * Phi**2 / 2) ** 2) / (2 * M)
    return {
        "families": rows(),
        "counterterm_occurrences": counterterms(),
        "scalar_refinements": sum(r.get("refinements", 0) for r in rows()),
        "fermionic_quartic_families": [r["id"] for r in rows()[4:]],
        "heavy_source_eliminated_action": eliminated,
        "source_dependent_action": s.expand(eliminated + G**2 * Phi**4 / (8 * M)),
        "checks": {
            "four_scalar_families": len(rows()[:4]) - 4,
            "all_192_scalar_refinements": sum(r.get("refinements", 0) for r in rows())
            - 192,
            "four_fermionic_quartic_rows": len(rows()[4:]) - 4,
            "only_two_new_matching_occurrences": sum(r["extra"] for r in counterterms())
            - 2,
            "heavy_source_no_quartic_variation": s.diff(
                eliminated, J, Phi, Phi, Phi, Phi
            ),
            "heavy_source_only_mass_and_constant": s.expand(
                eliminated
                + G**2 * Phi**4 / (8 * M)
                + J**2 / (2 * M)
                + J * G * Phi**2 / (2 * M)
            ),
            "local_contact_no_forward_second_coefficient": s.diff(
                s.Symbol("local_contact"), s.Symbol("nu"), 2
            ),
        },
        "scope": "All fixed-hybrid four-point families with their assigned proper and overall references. Extra=True means a new coordinate-matching contribution, not another raw graph. H elimination is exact; source dependence has degree at most two in Phi. Phi parity excludes internal light tree exchange.",
    }
