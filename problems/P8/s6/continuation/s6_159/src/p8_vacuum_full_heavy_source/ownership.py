"""All order-two H one-point occurrences in the fixed physical-Phi model."""


def rows():
    return [
        {
            "id": "scalar_OS_covariance_tadpole",
            "owner": "differentiated_scalar_sunset_and_OS_references",
            "status": "BOUNDED",
        },
        {
            "id": "fermionic_OS_covariance_tadpole",
            "owner": "S6.148",
            "status": "BOUNDED",
        },
        {
            "id": "proper_cubic_MS_counterterm_tadpole",
            "owner": "this_checkpoint",
            "status": "BOUNDED",
        },
        {
            "id": "first_Phi_coordinate_variation_of_regulated_J1",
            "owner": "this_checkpoint",
            "status": "BOUNDED",
        },
        {
            "id": "scalar_local_L_tadpole_and_mass_reference",
            "owner": "fixed_reference_OS_cancellation",
            "status": "CANCELLED",
        },
        {
            "id": "heavy_reducible_source_occurrences",
            "owner": "exact_Gaussian_stationary_identity",
            "status": "CANCELLED",
        },
        {
            "id": "heavy_MS_mass_and_fermion_MS_counterterm_vacuum_derivatives",
            "owner": "independent_of_H_at_Phi_zero",
            "status": "ZERO",
        },
        {
            "id": "second_mass_cubic_and_field_counterterms_inside_one_loop",
            "owner": "first_enter_at_order_three",
            "status": "HIGHER_ORDER",
        },
    ]


def validate_rows(value):
    if value != rows():
        raise ValueError("The complete heavy-source ownership ledger changed")
    return True


def data():
    return {
        "rows": rows(),
        "fixed_source_condition": "At H=Phi=0, the exact Gaussian H equation is J_B+G_B <Phi_B^2>/2=0. In canonical coefficient language this includes the separate proper cubic MS counterterm times the free light tadpole. It does not supply a new heavy-field normalization condition.",
        "checks": {
            "four_nonzero_source_classes": sum(
                row["status"] == "BOUNDED" for row in rows()
            )
            - 4,
            "two_assigned_cancellations": sum(
                row["status"] == "CANCELLED" for row in rows()
            )
            - 2,
            "eight_complete_source_ownership_classes": len(rows()) - 8,
        },
        "scope": "Full stationary H one-point reference at order two for the named Phi/vacuum observables. No global potential, external-fermion canonical dictionary, derivative-coordinate source dictionary or physical truncation is implied.",
    }
