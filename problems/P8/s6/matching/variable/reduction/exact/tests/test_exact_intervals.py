"""Independent arithmetic and derivative controls for the rational box proof."""

import json
from fractions import Fraction as Q

import pytest
from p8_exact_stationary import intervals as box


@pytest.mark.parametrize(("left", "right", "product"), [
    ((-2, 3), (-4, 5), (-12, 15)),
    ((-3, -2), (-5, -4), (8, 15)),
    ((2, 3), (4, 5), (8, 15)),
    ((0, 0), (-4, 5), (0, 0)),
])
def test_interval_product_extrema(left, right, product):
    assert box.Interval(*left) * box.Interval(*right) == box.Interval(*product)


def test_interval_reciprocal_and_integer_powers():
    assert box.Interval(-4, -2).inverse() == box.Interval(Q(-1, 2), Q(-1, 4))
    assert box.Interval(-2, 3)**2 == box.Interval(0, 9)
    assert box.Interval(-2, 3)**3 == box.Interval(-8, 27)
    assert box.Interval(2, 4)**(-2) == box.Interval(Q(1, 16), Q(1, 4))
    assert box.Interval(-2, 3)**0 == box.Interval(1)
    assert 1 - box.Interval(2, 3) == box.Interval(-2, -1)
    assert 2 / box.Interval(2, 4) == box.Interval(Q(1, 2), 1)


@pytest.mark.parametrize("value", [True, False, 0.0, float("nan"), float("inf"), "1/2", None])
def test_strict_exact_inputs(value):
    with pytest.raises(TypeError):
        box.rational(value)
    with pytest.raises(TypeError):
        box.Interval(value)
    with pytest.raises(TypeError):
        _ = box.Interval(1) + value


def test_invalid_interval_operations_are_rejected():
    with pytest.raises(ValueError):
        box.Interval(2, 1)
    with pytest.raises(ZeroDivisionError):
        box.Interval(-1, 1).inverse()
    with pytest.raises(ZeroDivisionError):
        _ = box.Interval(0, 1)**(-1)
    for exponent in (True, Q(1, 2), 2.0):
        with pytest.raises(TypeError):
            _ = box.Interval(1, 2)**exponent


def test_outward_rounding_handles_both_negative_endpoints():
    value = box.Interval(Q(-17, 100), Q(-13, 100))
    assert box.outward(value, 10) == box.Interval(Q(-1, 5), Q(-1, 10))
    assert box.outward(box.Interval(Q(-1, 10), Q(1, 10)), 10) == box.Interval(Q(-1, 10), Q(1, 10))
    for bad in (True, Q(10), 10.0):
        with pytest.raises(TypeError):
            box.outward(value, bad)
    for bad in (0, -1):
        with pytest.raises(ValueError):
            box.outward(value, bad)


def test_full_box_reproduces_exact_rational_headline_bounds():
    data = box.calibration()
    bounds = data["outward_enclosures"]
    assert bounds["T"] == box.Interval(Q(11_864_487, 10**6), Q(12_164_087, 10**6))
    assert bounds["T_zeta"] == box.Interval(Q(-5, 10**6), Q(1627, 10**6))
    assert bounds["N"] == box.Interval(Q(1_981_609, 10**6), Q(2_042_137, 10**6))
    assert all(value > 0 for value in data["strict_margins"].values())
    assert all(box.checks().values())


def test_cube_root_endpoints_are_proved_not_rounded():
    data = box.enclose()
    lower, upper = Q(199, 100), Q(201, 100)
    assert lower**3 < data["b_cubed"].lo
    assert data["b_cubed"].hi < upper**3
    assert lower > Q(19, 10) and upper < Q(21, 10)
    assert data["F"].hi < 0
    assert data["b_total_v"].hi < 0


@pytest.mark.parametrize(("key", "bad"), [
    ("v", box.Interval(Q(-1, 10**6), 0)),
    ("v", box.Interval(0, Q(2, 10_000))),
    ("delta", box.Interval(-1, 0)),
    ("delta", box.Interval(0, Q(2, 100))),
    ("zeta", box.Interval(10, 13)),
    ("zeta", box.Interval(11, 14)),
])
def test_parameter_domain_is_not_silently_enlarged(key, bad):
    domain = box.fixed_box()
    domain[key] = bad
    with pytest.raises(ValueError):
        box.enclose(domain)


def test_box_schema_and_exact_interval_entry_guards():
    with pytest.raises(TypeError):
        box.enclose(False)
    with pytest.raises(ValueError):
        box.enclose({"v": box.Interval(0)})
    domain = box.fixed_box()
    domain["c"] = box.Interval(2)
    with pytest.raises(ValueError):
        box.enclose(domain)
    domain = box.fixed_box()
    domain["v"] = (0, Q(1, 10_000))
    with pytest.raises(TypeError):
        box.enclose(domain)


class Jet:
    """Separate exact first-order dual arithmetic for point derivative tests."""

    def __init__(self, value, derivative=0):
        self.value, self.derivative = Q(value), Q(derivative)

    def __add__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value + other.value, self.derivative + other.derivative)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, -self.derivative)

    def __sub__(self, other):
        return self + -other if isinstance(other, Jet) else self + Jet(-other)

    def __rsub__(self, other):
        return other + -self

    def __mul__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value * other.value,
                   self.derivative * other.value + self.value * other.derivative)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value / other.value,
                   (self.derivative * other.value - self.value * other.derivative) / other.value**2)

    def __rtruediv__(self, other):
        return Jet(other) / self

    def __pow__(self, exponent):
        return Jet(self.value**exponent, exponent*self.value**(exponent - 1)*self.derivative)


@pytest.mark.parametrize("c", [Q(2), Q(2001, 1000), Q(201, 100)])
def test_partials_against_independent_exact_forward_jets_at_center(c):
    v, d = Jet(0, 1), Jet(1, 1)
    zeta = 3*(c + 2)**2/(2*c)
    J = c*d**4*(1 - 7*v) + 12*v
    Jv = c*d**3*(-3 - 35*v) + 12
    Q0 = 32*(1 - v)/(c*d**14)
    R = d**8*J/(8*c*(1 - v))
    L = 8/d + Jv/J + 1/(1 - v)
    D = c - 2/d**4
    H = D*(Jv/J - 6/d) - 8/d**5
    F = 2*(v*zeta*H - L)/3
    b3 = (1 - v*D*zeta)/R
    assert b3.value == 8
    b = Jet(2, b3.derivative/12)
    T = 3*b*F**2/(2*Q0)
    exact = {"J": J.value, "J_v": J.derivative, "J_vv": Jv.derivative,
             "Q": Q0.value, "Q_v_over_Q": Q0.derivative/Q0.value,
             "R": R.value, "L": L.value, "L_v": L.derivative,
             "D": D.value, "D_v": D.derivative, "H": H.value,
             "H_v": H.derivative, "F": F.value, "F_v": F.derivative,
             "b_v": b.derivative, "T": T.value, "T_v": T.derivative}
    actual = box.enclose({"v": box.Interval(0), "delta": box.Interval(c - 2),
                          "zeta": box.Interval(zeta)})
    assert T.value == zeta
    assert b.derivative == -c*(c + 2)
    assert 2*b.derivative/F.value == c**2/2
    for name, value in exact.items():
        assert actual[name].contains(value), name


def test_implicit_lapse_correction_cannot_be_omitted():
    data = box.enclose({"v": box.Interval(Q(1, 10_000)),
                        "delta": box.Interval(Q(1, 100)), "zeta": box.Interval(12)})
    assert data["b_zeta"].hi < 0
    assert data["T_v"].lo > 0
    assert data["F"].hi < 0
    assert data["lapse_implicit_correction"].lo > 0


def test_report_is_exact_serializable_and_keeps_boundary_scope():
    report = box.report()
    json.dumps(report, sort_keys=True, allow_nan=False)

    def no_floats(value):
        if isinstance(value, dict):
            return all(no_floats(item) for item in value.values())
        if isinstance(value, list):
            return all(no_floats(item) for item in value)
        return not isinstance(value, float)

    assert no_floats(report)
    assert report["literal_delta_positive"]
    assert report["delta_zero_is_only_desingularized_limit"]
    assert not report["full_parent_solution"]
    assert not report["tensor_Green_inverse_selected"]
    assert report["u_half_width"] == "1/100"
